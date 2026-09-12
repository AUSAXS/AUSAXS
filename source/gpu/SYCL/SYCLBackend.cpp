// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/* The kernel is structured around three decisions:

 - Bins are 64-bit fixed-point integers, because every contribution reaches its bin through an atomic add and only the integer one is a 
   machine instruction. f32 is imprecise, and f64 has only little support on modern consumer GPUs. 

 - Contributions go into a workgroup-local histogram that is flushed once per tile. Only the first local_bins bins fit in workgroup memory; 
   the tail falls back to global memory, which is rare in practice since it covers the first local_bins*bin_width angstrom.

 - Work is split into a uniform tile grid over both index ranges, and each thread keeps rows_per_thread atoms of the first range in registers 
   while streaming the second. The streamed atom is the same for every thread in the workgroup, so it is read once into scalar registers and 
   its cost is amortized over all the atoms held in registers. This is intended to minimize memory traffic. 
*/

#include <gpu/GPUBackendABI.h>

#include <sycl/sycl.hpp>

#include <array>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace ausaxs::gpu::abi;

namespace {
    using i64 = std::int64_t;

    constexpr int workgroup_size = 256;
    constexpr int rows_per_thread = 4;                          // atoms of the first range held in registers
    constexpr int tile_size = workgroup_size*rows_per_thread;   // atoms covered by one tile in either direction

    /**
     * @brief How many leading bins to keep in workgroup memory. This must be strictly less than the local memory limit. 
     */
    int local_bin_capacity(std::size_t local_mem, int bin_count, bool weighted) {
        const std::size_t per_bin = (weighted ? 3 : 1)*sizeof(i64);
        return static_cast<int>(sycl::min(static_cast<std::size_t>(bin_count), local_mem/per_bin));
    }

    /**
     * @brief How many output histograms to make room for before any are asked for.
     */
    constexpr int default_slots = 256;

    // chosen so that a bin's fractional part keeps 16 bits of precision while the whole accumulated range stays far below the int64 overflow point.
    constexpr float fixed_scale = 65536.f;

    // Convert a value to its fixed-point representation.
    inline i64 to_fixed(float value) {
        return static_cast<i64>(sycl::rint(value*fixed_scale));
    }

    // Convert a fixed-point accumulator back to the real value the ABI reports.
    inline double from_fixed(i64 value) {
        return static_cast<double>(value)/fixed_scale;
    }

    template<sycl::access::address_space space, sycl::memory_scope scope>
    inline void atomic_add(i64& target, i64 value) {
        sycl::atomic_ref<i64, sycl::memory_order::relaxed, scope, space> ref{target};
        ref.fetch_add(value);
    }

    inline void add_local(i64& target, i64 value) {
        atomic_add<sycl::access::address_space::local_space, sycl::memory_scope::work_group>(target, value);
    }

    inline void add_global(i64& target, i64 value) {
        atomic_add<sycl::access::address_space::global_space, sycl::memory_scope::device>(target, value);
    }

    // A queued calculation as the kernel sees it.
    struct DeviceJob {
        const sycl::float4* a1;
        const sycl::float4* a2; // equal to a1 for self-correlations
        std::uint32_t n1, n2;
        std::uint32_t scaling;
        std::uint32_t slot;
    };

    // One tile of one job, i.e. the work of a single workgroup.
    struct TileRef {
        std::uint32_t job;
        std::uint32_t i_block;
        std::uint32_t j_block;
    };

    /**
     * @brief Evaluate one tile of the distance matrix into the workgroup-local histogram.
     *
     * @tparam weighted Whether to track the bin count and the summed distance alongside the value.
     * @tparam diagonal Whether the tile lies on the diagonal of a self-correlation, in which case only the pairs above the diagonal are evaluated. 
     */
    template<bool weighted, bool diagonal>
    struct Kernel {
        const TileRef* tiles;
        const DeviceJob* jobs;
        i64* histograms;         // n_slots*bin_count values, or the same number of {value, count, center} triples
        float inv_width;
        int bin_count;
        int local_bins; // leading bins held in workgroup memory; never exceeds bin_count

        // workgroup-local accumulators; count and center are unused when not weighted
        using local_array = sycl::local_accessor<i64, 1>;
        local_array local_value, local_count, local_center;

        /**
         * @brief Accumulate one contribution, in workgroup memory if the bin is held there.
         */
        void add_to_bin(int bin, float value, float center, std::uint32_t count, i64* global) const {
            if (bin < local_bins) {
                add_local(local_value[bin], to_fixed(value));
                if constexpr (weighted) {
                    add_local(local_count[bin], count);
                    add_local(local_center[bin], to_fixed(center));
                }
            } else {
                if constexpr (weighted) {
                    add_global(global[3*bin + 0], to_fixed(value));
                    add_global(global[3*bin + 1], count);
                    add_global(global[3*bin + 2], to_fixed(center));
                } else {
                    add_global(global[bin], to_fixed(value));
                }
            }
        }

        void operator()(sycl::nd_item<1> item) const {
            const int lid = static_cast<int>(item.get_local_id(0));
            const auto group = item.get_group();

            for (int bin = lid; bin < local_bins; bin += workgroup_size) {
                local_value[bin] = 0;
                if constexpr (weighted) {
                    local_count[bin] = 0;
                    local_center[bin] = 0;
                }
            }
            sycl::group_barrier(group);

            const TileRef tile = tiles[item.get_group(0)];
            const DeviceJob job = jobs[tile.job];
            i64* global = histograms + static_cast<std::size_t>(job.slot)*bin_count*(weighted ? 3 : 1);

            // the atoms this thread keeps in registers for the whole tile
            std::array<sycl::float4, rows_per_thread> held;
            std::array<int, rows_per_thread> index;
            std::array<bool, rows_per_thread> active;
            for (int k = 0; k < rows_per_thread; ++k) {
                const int local = lid + k*workgroup_size;
                index[k] = static_cast<int>(tile.i_block)*tile_size + local;
                active[k] = index[k] < static_cast<int>(job.n1);
                held[k] = active[k] ? job.a1[index[k]] : sycl::float4{0.f, 0.f, 0.f, 0.f};
            }

            // every contribution counts the pair twice, since each unordered pair is visited once
            const float scale = 2.f*static_cast<float>(job.scaling);
            const std::uint32_t count = 2*job.scaling;

            const int j_start = static_cast<int>(tile.j_block)*tile_size;
            const int j_end = sycl::min(j_start + tile_size, static_cast<int>(job.n2));
            for (int j = j_start; j < j_end; ++j) {
                const sycl::float4 other = job.a2[j]; // identical across the workgroup, so read into scalar registers
                const float w = scale*other.w();
                for (int k = 0; k < rows_per_thread; ++k) {
                    if constexpr (diagonal) {
                        if (!active[k] || j <= index[k]) {continue;}
                    } else {
                        if (!active[k]) {continue;}
                    }
                    const float dx = held[k].x() - other.x();
                    const float dy = held[k].y() - other.y();
                    const float dz = held[k].z() - other.z();
                    const float d = sycl::sqrt(dx*dx + dy*dy + dz*dz);
                    add_to_bin(
                        static_cast<int>(sycl::rint(inv_width*d)),
                        w*held[k].w(), scale*d, count, global
                    );
                }
            }

            sycl::group_barrier(group);
            for (int bin = lid; bin < local_bins; bin += workgroup_size) {
                if constexpr (weighted) {
                    if (i64 v = local_value[bin]; v != 0) {add_global(global[3*bin + 0], v);}
                    if (i64 c = local_count[bin]; c != 0) {add_global(global[3*bin + 1], c);}
                    if (i64 c = local_center[bin]; c != 0) {add_global(global[3*bin + 2], c);}
                } else {
                    if (i64 v = local_value[bin]; v != 0) {add_global(global[bin], v);}
                }
            }
        }
    };

    /**
     * @brief Bump allocator for the buffers one calculation submits.
     *
     * A submit may not wait for the work it queues, so the buffers it fills stay in use for an unknown time afterwards and cannot be recycled 
     * by the next one. Allocations therefore only move forward within a calculation, and the whole arena is rewound in begin(), by which point
     * finish() has waited for everything that could still be reading it. The slabs themselves are kept and reused, so a steady-state run stops 
     * allocating entirely.
     */
    class Arena {
        public:
            void* allocate(sycl::queue& queue, std::size_t bytes) {
                if (bytes == 0) {return nullptr;}
                bytes = (bytes + 255)/256*256; // keep every allocation comfortably aligned

                for (; slab < slabs.size(); ++slab, used = 0) {
                    if (used + bytes <= slabs[slab].size) {break;}
                }
                if (slab == slabs.size()) {
                    const std::size_t size = sycl::max(bytes, default_slab);
                    void* pointer = sycl::malloc_device(size, queue);
                    if (pointer == nullptr) {throw std::runtime_error("sycl_backend: out of device memory");}
                    slabs.push_back(Slab{.ptr=pointer, .size=size});
                    used = 0;
                }

                void* result = static_cast<char*>(slabs[slab].ptr) + used;
                used += bytes;
                return result;
            }

            // Make the whole arena available again. Only safe once the device is idle.
            void rewind() {
                slab = 0;
                used = 0;
            }

            void release(sycl::queue& queue) {
                for (auto& s : slabs) {sycl::free(s.ptr, queue);}
                slabs.clear();
                rewind();
            }

        private:
            struct Slab {void* ptr; std::size_t size;};
            static constexpr std::size_t default_slab = 8*1024*1024;

            std::vector<Slab> slabs;
            std::size_t slab = 0, used = 0;
    };

    /**
     * @brief The device and the buffers shared by all calculations.
     *
     * Acquiring a device and compiling the kernels is expensive, so both are done once and kept for the lifetime of the process.
     */
    class Context {
        public:
            static Context* get() {
                static std::unique_ptr<Context> instance = [] () -> std::unique_ptr<Context> {
                    try {
                        return std::unique_ptr<Context>(new Context());
                    } catch (...) {
                        return nullptr;
                    }
                }();
                return instance.get();
            }

            sycl::queue queue;
            std::string name;
            Arena arena;
            std::size_t local_mem = 0; // workgroup memory the device offers a single workgroup

            // state of the calculation currently being built, reset by begin()
            bool active = false;
            bool weighted = false;
            int bin_count = 0;
            int local_bins = 0; // derived from local_mem and bin_count; see local_bin_capacity
            float inv_width = 0;
            std::unordered_map<const float*, const sycl::float4*> uploaded; // host -> device, deduplicated
            std::vector<sycl::event> pending;                               // uploads the next kernel must wait for
            sycl::event zeroing;                                            // begin_on clearing the slots the last calculation used

            /**
             * @brief One past the highest slot any submit has accumulated into since it was cleared.
             */
            int dirty_slots = 0;

            /**
             * @brief The output histograms, slot-major.
             *
             * Kept across calculations and only ever grown, so the reallocation below happens a
             * handful of times in a process rather than once per calculation.
             */
            i64* histograms = nullptr;
            int slot_capacity = 0;

            /**
             * @brief Make room for @a slots histograms, preserving what is already accumulated.
             *
             * The caller may name a slot it has never used before in any submit, so this can happen with kernels in flight — hence the wait. 
             * Growth is slot-major and therefore a plain prefix copy. It is the only synchronisation between begin() and finish().
             */
            void ensure_slots(int slots) {
                if (slots <= slot_capacity) {return;}
                slots = sycl::max(slots, 2*slot_capacity);

                const std::size_t stride = weighted ? 3 : 1;
                const std::size_t values = static_cast<std::size_t>(slots)*bin_count*stride;
                auto* fresh = static_cast<i64*>(sycl::malloc_device(values*sizeof(i64), queue));
                if (fresh == nullptr) {throw std::runtime_error("sycl_backend: out of device memory");}
                queue.fill(fresh, i64{0}, values).wait();

                if (histograms != nullptr) {
                    queue.wait(); // in-flight kernels are still accumulating into the old buffer
                    const std::size_t old = static_cast<std::size_t>(slot_capacity)*bin_count*stride;
                    queue.memcpy(fresh, histograms, old*sizeof(i64)).wait();
                    sycl::free(histograms, queue);
                }
                histograms = fresh;
                slot_capacity = slots;
                zeroing = {}; // everything above was waited for, and the buffer it referred to is gone
            }

        private:
            Context() : queue(sycl::gpu_selector_v) {
                name = queue.get_device().get_info<sycl::info::device::name>();
                local_mem = queue.get_device().get_info<sycl::info::device::local_mem_size>();
                warmup();
            }

            // Compile all kernel variants on a trivial input, so the first real calculation does not pay for it.
            void warmup();
    };

    /**
     * @brief Upload any coordinate set this calculation has not seen yet, and return the device pointers in job order.
     *
     * The deduplication map lives for the whole calculation rather than one submit, so a coordinate set shared by jobs in different submits 
     * is still uploaded only once - which is the common case for the symmetry managers, where one body's coordinates appear in many jobs. 
     * The uploads are asynchronous; their events are collected so the kernels of this submit can wait on them without the host having to.
     */
    std::vector<std::pair<const sycl::float4*, const sycl::float4*>> upload_coordinates(
        Context& context, const Job* jobs, int n_jobs
    ) {
        auto upload = [&] (const float* host, std::uint32_t n) -> const sycl::float4* {
            if (auto it = context.uploaded.find(host); it != context.uploaded.end()) {return it->second;}
            const std::size_t bytes = static_cast<std::size_t>(n)*sizeof(sycl::float4);
            auto* device = static_cast<sycl::float4*>(context.arena.allocate(context.queue, bytes));
            context.pending.push_back(context.queue.memcpy(device, host, bytes));
            context.uploaded[host] = device;
            return device;
        };

        std::vector<std::pair<const sycl::float4*, const sycl::float4*>> result(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const sycl::float4* a1 = upload(jobs[i].a1, jobs[i].n1);
            result[i] = {a1, (jobs[i].a2 != nullptr) ? upload(jobs[i].a2, jobs[i].n2) : a1};
        }
        return result;
    }

    /**
     * @brief Split the queued jobs into the tiles that make up one workgroup of work each.
     *
     * Self-correlations only enumerate the tiles on and above the diagonal, which are returned separately as they need the extra predicate.
     */
    void build_tiles(
        const Job* jobs, int n_jobs, std::vector<TileRef>& regular, std::vector<TileRef>& diagonal
    ) {
        auto blocks = [] (std::uint32_t n) {return (n + tile_size - 1)/tile_size;};
        for (int i = 0; i < n_jobs; ++i) {
            const auto job = static_cast<std::uint32_t>(i);
            if (jobs[i].a2 == nullptr) {
                const std::uint32_t n = blocks(jobs[i].n1);
                for (std::uint32_t p = 0; p < n; ++p) {
                    diagonal.emplace_back(TileRef{.job=job, .i_block=p, .j_block=p});
                    for (std::uint32_t q = p + 1; q < n; ++q) {
                        regular.emplace_back(TileRef{.job=job, .i_block=p, .j_block=q});
                    }
                }
            } else {
                for (std::uint32_t p = 0; p < blocks(jobs[i].n1); ++p) {
                    for (std::uint32_t q = 0; q < blocks(jobs[i].n2); ++q) {
                        regular.emplace_back(TileRef{.job=job, .i_block=p, .j_block=q});
                    }
                }
            }
        }
    }

    template<bool weighted, bool diagonal>
    void launch(
        Context& context, const TileRef* tiles, std::size_t n_tiles, const DeviceJob* jobs,
        i64* histograms, float inv_width, int bin_count, int local_bins, // NOLINT
        const std::vector<sycl::event>& wait_for
    ) {
        if (n_tiles == 0) {return;}
        context.queue.submit([&] (sycl::handler& handler) {
            handler.depends_on(wait_for);
            handler.depends_on(context.zeroing);

            // the two the unweighted kernel never touches are left empty rather than given a token
            // element, so that the allocation is exactly what local_bin_capacity budgeted for
            const auto bins = static_cast<std::size_t>(local_bins);
            const auto shared = sycl::range<1>(weighted ? bins : 0);
            sycl::local_accessor<i64, 1> value{sycl::range<1>(bins), handler};
            sycl::local_accessor<i64, 1> count{shared, handler};
            sycl::local_accessor<i64, 1> center{shared, handler};
            handler.parallel_for(
                sycl::nd_range<1>(n_tiles*workgroup_size, workgroup_size),
                Kernel<weighted, diagonal>{
                    tiles, jobs, histograms, inv_width, bin_count,
                    local_bins, value, count, center
                }
            );
        });
    }

    // Queue one batch of correlations. Does not wait for them; see the ABI header.
    template<bool weighted>
    void submit_on(Context& context, const Job* jobs, int n_jobs) {
        if (n_jobs == 0) {return;}

        std::uint32_t max_slot = 0;
        for (int i = 0; i < n_jobs; ++i) {max_slot = sycl::max(max_slot, jobs[i].slot);}
        context.ensure_slots(static_cast<int>(max_slot) + 1);
        context.dirty_slots = sycl::max(context.dirty_slots, static_cast<int>(max_slot) + 1);

        const auto coordinates = upload_coordinates(context, jobs, n_jobs);
        std::vector<DeviceJob> device_jobs(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const bool self = jobs[i].a2 == nullptr;
            device_jobs[i] = DeviceJob{
                .a1=coordinates[i].first, .a2=coordinates[i].second, .n1=jobs[i].n1,
                .n2=self ? jobs[i].n1 : jobs[i].n2, .scaling=jobs[i].scaling, .slot=jobs[i].slot
            };
        }

        std::vector<TileRef> regular, diagonal;
        build_tiles(jobs, n_jobs, regular, diagonal);

        // every buffer below is arena memory belonging to this submit alone, so the kernels may
        // still be reading it when the next submit arrives
        const std::size_t job_bytes = device_jobs.size()*sizeof(DeviceJob);
        const std::size_t tile_bytes = (regular.size() + diagonal.size())*sizeof(TileRef);
        auto* device_job_array = static_cast<DeviceJob*>(context.arena.allocate(context.queue, job_bytes));
        auto* device_tiles = static_cast<TileRef*>(context.arena.allocate(context.queue, tile_bytes));

        context.pending.push_back(context.queue.memcpy(device_job_array, device_jobs.data(), job_bytes));
        if (!regular.empty()) {
            context.pending.push_back(
                context.queue.memcpy(device_tiles, regular.data(), regular.size()*sizeof(TileRef))
            );
        }
        if (!diagonal.empty()) {
            context.pending.push_back(context.queue.memcpy(
                device_tiles + regular.size(), diagonal.data(), diagonal.size()*sizeof(TileRef)
            ));
        }

        launch<weighted, false>(
            context, device_tiles, regular.size(), device_job_array,
            context.histograms, context.inv_width, context.bin_count, context.local_bins, context.pending
        );
        launch<weighted, true>(
            context, device_tiles + regular.size(), diagonal.size(), device_job_array,
            context.histograms, context.inv_width, context.bin_count, context.local_bins, context.pending
        );

        // the uploads are now sequenced before the kernels that read them; a later submit brings its
        // own, and must not make its kernels wait on this one's as well
        context.pending.clear();
    }

    /**
     * @brief Wait for everything submitted since begin(), then copy the raw fixed-point histograms
     *        into a host staging buffer so they can be converted; shared by both finish variants
     */
    std::vector<i64> read_raw(Context& context, int n_slots, int stride) {
        // check_finish has already bounded n_slots by what was submitted, so nothing has to be allocated or cleared here to answer it
        context.queue.wait_and_throw();

        std::vector<i64> raw(static_cast<std::size_t>(n_slots)*context.bin_count*stride);
        if (!raw.empty()) {context.queue.memcpy(raw.data(), context.histograms, raw.size()*sizeof(i64)).wait();}
        context.active = false;
        return raw;
    }

    // Wait for everything submitted since begin() and convert the histograms into @a out.
    void finish_unweighted_on(Context& context, int n_slots, double* out) {
        const auto raw = read_raw(context, n_slots, 1);
        for (std::size_t i = 0; i < raw.size(); ++i) {out[i] = from_fixed(raw[i]);}
    }

    // Wait for everything submitted since begin() and convert the histograms into @a out.
    void finish_weighted_on(Context& context, int n_slots, WeightedBin* out) {
        const auto raw = read_raw(context, n_slots, 3);
        for (std::size_t i = 0; i < raw.size()/3; ++i) {
            out[i].value  = from_fixed(raw[3*i + 0]);
            out[i].count  = raw[3*i + 1];
            out[i].center = from_fixed(raw[3*i + 2]);
        }
    }

    //  Start a calculation on an idle device, discarding whatever the last one left behind.
    void begin_on(Context& context, int bin_count, float inv_width, bool weighted) {
        context.queue.wait_and_throw(); // nothing may still be reading the arena we are about to rewind
        context.arena.rewind();
        context.uploaded.clear();
        context.pending.clear();

        // the layout of the output depends on both, so a change invalidates what is already allocated
        if (context.bin_count != bin_count || context.weighted != weighted) {
            if (context.histograms != nullptr) {sycl::free(context.histograms, context.queue);}
            context.histograms = nullptr;
            context.slot_capacity = 0;
            context.dirty_slots = 0; // the replacement ensure_slots allocates below arrives cleared
        }
        context.bin_count = bin_count;
        context.local_bins = local_bin_capacity(context.local_mem, bin_count, weighted);
        context.inv_width = inv_width;
        context.weighted = weighted;
        context.active = true;

        context.ensure_slots(sycl::max(context.slot_capacity, default_slots));

        // only what the previous calculation accumulated into needs clearing; every slot above that is untouched since ensure_slots allocated 
        // it, and so already zero. the fill is not waited for - the host has coordinates to upload and tiles to build first - so instead every
        // kernel of this calculation is made to depend on it.
        const std::size_t stride = weighted ? 3 : 1;
        const std::size_t values = static_cast<std::size_t>(context.dirty_slots)*bin_count*stride;
        context.zeroing = values == 0 ? sycl::event{} : context.queue.fill(context.histograms, i64{0}, values);
        context.dirty_slots = 0;
    }

    void Context::warmup() {
        // one self and one cross job, so both the diagonal and the regular kernel are compiled
        const std::vector<float> atoms(4*tile_size, 1.f);
        const std::array<Job, 2> jobs = {
            Job{.a1=atoms.data(), .a2=nullptr, .n1=tile_size, .n2=0, .scaling=1, .slot=0},
            Job{.a1=atoms.data(), .a2=atoms.data(), .n1=tile_size, .n2=tile_size, .scaling=1, .slot=0}
        };
        begin_on(*this, 8, 1.f, false);
        submit_on<false>(*this, jobs.data(), 2);
        std::vector<double> out_unweighted(8);
        finish_unweighted_on(*this, 1, out_unweighted.data());

        begin_on(*this, 8, 1.f, true);
        submit_on<true>(*this, jobs.data(), 2);
        std::vector<WeightedBin> out_weighted(8);
        finish_weighted_on(*this, 1, out_weighted.data());
    }

}

namespace {
    // the message belonging to the most recent failure on this thread, as returned by last_error()
    thread_local std::string error_message;

    /**
     * @brief Translate any failure into a status code. This is required to avoid exceptions propagating across the ABI boundary. 
     */
    template<typename work_t>
    Status guarded(work_t&& work) {
        error_message.clear();
        try {
            Context* context = Context::get();
            if (!context) {
                error_message = "sycl_backend: no usable device";
                return Status::no_device;
            }
            work(*context);
            return Status::ok;
        } catch (const sycl::exception& e) {
            error_message = std::string("sycl_backend: ") + e.what();
            return e.code() == sycl::errc::memory_allocation ? Status::out_of_memory : Status::kernel_failed;
        } catch (const std::bad_alloc&) {
            error_message = "sycl_backend: out of host memory";
            return Status::out_of_memory;
        } catch (const std::exception& e) {
            error_message = std::string("sycl_backend: ") + e.what();
            // the only runtime_errors thrown here are allocation failures
            return Status::out_of_memory;
        } catch (...) {
            error_message = "sycl_backend: unknown failure";
            return Status::kernel_failed;
        }
    }

    // Reject a call that cannot be served, before it reaches the device.
    Status check_finish(Context* context, bool weighted, int n_slots) {
        if ((context == nullptr) || !context->active) {
            error_message = "sycl_backend: finish without a matching begin";
            return Status::invalid_input;
        }
        if (context->weighted != weighted) {
            error_message = "sycl_backend: finish does not match the weighted flag given to begin";
            return Status::invalid_input;
        }
        if (n_slots < 0) {
            error_message = "sycl_backend: negative slot count";
            return Status::invalid_input;
        }
        // nothing was reserved for a slot no job named, so there is no result to hand back for one
        if (context->dirty_slots < n_slots) {
            error_message =
                "sycl_backend: finish asks for " + std::to_string(n_slots) + " slots, but only " +
                std::to_string(context->dirty_slots) + " were submitted for";
            return Status::invalid_input;
        }
        return Status::ok;
    }
}

std::int32_t ausaxs::gpu::abi::ausaxs_gpu_abi_version() {
    return ausaxs::gpu::abi::version;
}

bool ausaxs::gpu::abi::ausaxs_gpu_available() {
    return Context::get() != nullptr;
}

const char* ausaxs::gpu::abi::ausaxs_gpu_device_name() {
    Context* context = Context::get();
    return (context != nullptr) ? context->name.c_str() : "none";
}

const char* ausaxs::gpu::abi::ausaxs_gpu_last_error() {
    return error_message.c_str();
}

Status ausaxs::gpu::abi::ausaxs_gpu_begin(std::int32_t bin_count, float inv_width, bool weighted) {
    if (bin_count <= 0) {
        error_message = "sycl_backend: bin count must be positive";
        return Status::invalid_input;
    }
    return guarded([&] (Context& context) {begin_on(context, bin_count, inv_width, weighted);});
}

Status ausaxs::gpu::abi::ausaxs_gpu_submit(const Job* jobs, std::int32_t n_jobs) {
    if (n_jobs < 0) {
        error_message = "sycl_backend: negative job count";
        return Status::invalid_input;
    }
    for (std::int32_t i = 0; i < n_jobs; ++i) {
        if (jobs[i].a1 == nullptr) {
            error_message = "sycl_backend: job " + std::to_string(i) + " has no coordinates";
            return Status::invalid_input;
        }
    }
    if (Context* context = Context::get(); (context != nullptr) && !context->active) {
        error_message = "sycl_backend: submit without a matching begin";
        return Status::invalid_input;
    }
    return guarded([&] (Context& context) {
        if (context.weighted) {submit_on<true>(context, jobs, n_jobs);}
        else                  {submit_on<false>(context, jobs, n_jobs);}
    });
}

Status ausaxs::gpu::abi::ausaxs_gpu_finish_unweighted(std::int32_t n_slots, double* out) {
    if (auto status = check_finish(Context::get(), false, n_slots); status != Status::ok) {return status;}
    return guarded([&] (Context& context) {finish_unweighted_on(context, n_slots, out);});
}

Status ausaxs::gpu::abi::ausaxs_gpu_finish_weighted(std::int32_t n_slots, WeightedBin* out) {
    if (auto status = check_finish(Context::get(), true, n_slots); status != Status::ok) {return status;}
    return guarded([&] (Context& context) {finish_weighted_on(context, n_slots, out);});
}
