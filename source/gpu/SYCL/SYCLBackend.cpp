// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

/* The kernel is structured around three decisions:

 - Each pair is accumulated as a 64-bit fixed-point integer in workgroup memory, and each tile's sums are converted to f64 and added to the
   global histogram once, when the tile is flushed. Native integer atomics are the cheapest add on every device, while per-pair f64 work is
   slow wherever f64 is not a first-class type, which includes most consumer GPUs. The fixed-point scale only has to keep one tile's sums
   from overflowing, so it is chosen per job from the largest possible term and is finer than the precision of the float inputs. The job's
   pair factor is applied at the flush, so it never enters that bound, and the global f64 sums have no practical limit.

 - Contributions go into a workgroup-local histogram that is flushed once per tile. Only the first local_bins bins fit in workgroup memory; 
   the tail falls back to global memory, which is rare in practice since it covers the first local_bins*bin_width angstrom.

 - Work is split into a uniform tile grid over both index ranges, and each thread keeps rows_per_thread atoms of the first range in registers 
   while streaming the second. The streamed atom is the same for every thread in the workgroup, so it is read once into scalar registers and 
   its cost is amortized over all the atoms held in registers. This is intended to minimize memory traffic. 
*/

#include <gpu/GPUBackendABI.h>

#include <sycl/sycl.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace ausaxs::gpu::abi;

namespace {
    using i64 = std::int64_t;

    constexpr int workgroup_size = 256;
    constexpr int rows_per_thread = 4;                          // atoms of the first range held in registers
    constexpr int tile_size = workgroup_size*rows_per_thread;   // atoms covered by one tile in either direction

    /**
     * @brief One output bin, laid out exactly as the ABI returns it, so the finished histograms are copied back without conversion.
     */
    template<bool weighted>
    using bin_t = std::conditional_t<weighted, WeightedBin, double>;

    std::size_t bin_bytes(bool weighted) {
        return weighted ? sizeof(WeightedBin) : sizeof(double);
    }

    /**
     * @brief How many leading bins to keep in workgroup memory, each holding one fixed-point integer per tracked quantity.
     */
    int local_bin_capacity(std::size_t local_mem, int bin_count, bool weighted) {
        const std::size_t per_bin = (weighted ? 3 : 1)*sizeof(i64);
        return static_cast<int>(sycl::min(static_cast<std::size_t>(bin_count), local_mem/per_bin));
    }

    /**
     * @brief The fixed-point scale for terms no larger than @a bound: the largest power of two at which the sum of all tile_size^2 pairs of a
     *        tile still fits in an int64, even if every one of them lands in the same bin.
     */
    float fixed_scale(double bound) {
        if (bound <= 0) {return 1.f;} // every term is zero
        int exponent = 0;
        std::frexp(static_cast<double>(tile_size)*tile_size*bound, &exponent); // the tile sum is below 2^exponent
        return std::ldexp(1.f, sycl::min(62 - exponent, 100));
    }

    inline i64 to_fixed(float scaled) {
        return static_cast<i64>(sycl::rint(scaled));
    }

    /**
     * @brief How many output histograms to make room for before any are asked for.
     */
    constexpr int default_slots = 256;

    template<sycl::access::address_space space, sycl::memory_scope scope, typename T>
    inline void atomic_add(T& target, T value) {
        sycl::atomic_ref<T, sycl::memory_order::relaxed, scope, space> ref{target};
        ref.fetch_add(value);
    }

    template<typename T>
    inline void add_local(T& target, T value) {
        atomic_add<sycl::access::address_space::local_space, sycl::memory_scope::work_group>(target, value);
    }

    template<typename T>
    inline void add_global(T& target, T value) {
        atomic_add<sycl::access::address_space::global_space, sycl::memory_scope::device>(target, value);
    }

    // A queued calculation as the kernel sees it.
    struct DeviceJob {
        const sycl::float4* a1;
        const sycl::float4* a2;     // equal to a1 for self-correlations
        std::uint32_t n1, n2;
        std::uint32_t scaling;      // the caller's pair factor, applied as-is when a tile is flushed; see Job::scaling
        std::uint32_t slot;
        float value_scale;          // fixed-point scales of the workgroup-local sums; see fixed_scale
        float center_scale;
        double value_factor;        // scaling/value_scale, converting a flushed value sum to the real one
        double center_factor;       // scaling/center_scale
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
        bin_t<weighted>* histograms; // n_slots*bin_count bins, slot-major
        float inv_width;
        int bin_count;
        int local_bins; // leading bins held in workgroup memory; never exceeds bin_count

        // workgroup-local fixed-point sums of the unscaled terms; count and center are unused when not weighted
        sycl::local_accessor<i64, 1> local_value, local_count, local_center;

        /**
         * @brief Accumulate the pair at distance @a d with weight product @a value, in workgroup memory if its bin is held there.
         *        The rare pair beyond the local bins goes straight to the global histogram, already scaled.
         */
        void add_to_bin(int bin, float value, float d, const DeviceJob& job, bin_t<weighted>* global) const {
            if (bin < local_bins) {
                add_local(local_value[bin], to_fixed(value*job.value_scale));
                if constexpr (weighted) {
                    add_local(local_count[bin], i64{1});
                    add_local(local_center[bin], to_fixed(d*job.center_scale));
                }
            } else {
                const auto scaling = static_cast<double>(job.scaling);
                if constexpr (weighted) {
                    add_global(global[bin].value, scaling*value);
                    add_global(global[bin].count, static_cast<i64>(job.scaling));
                    add_global(global[bin].center, scaling*d);
                } else {
                    add_global(global[bin], scaling*value);
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
            bin_t<weighted>* global = histograms + static_cast<std::size_t>(job.slot)*bin_count;

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

            const int j_start = static_cast<int>(tile.j_block)*tile_size;
            const int j_end = sycl::min(j_start + tile_size, static_cast<int>(job.n2));
            for (int j = j_start; j < j_end; ++j) {
                const sycl::float4 other = job.a2[j]; // identical across the workgroup, so read into scalar registers
                const float w = other.w();
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
                    add_to_bin(static_cast<int>(sycl::rint(inv_width*d)), w*held[k].w(), d, job, global);
                }
            }

            sycl::group_barrier(group);
            for (int bin = lid; bin < local_bins; bin += workgroup_size) {
                if constexpr (weighted) {
                    if (i64 c = local_count[bin]; c != 0) {
                        add_global(global[bin].value, static_cast<double>(local_value[bin])*job.value_factor);
                        add_global(global[bin].count, c*static_cast<i64>(job.scaling));
                        add_global(global[bin].center, static_cast<double>(local_center[bin])*job.center_factor);
                    }
                } else {
                    if (i64 v = local_value[bin]; v != 0) {add_global(global[bin], static_cast<double>(v)*job.value_factor);}
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

    // A coordinate set on the device, with the largest weight magnitude it holds.
    struct Upload {
        const sycl::float4* device;
        float max_weight;
    };

    // A device the user may choose, with the identifier the ABI reports for it.
    struct Candidate {
        sycl::device device;
        std::string id;
    };

    /**
     * @brief Name of the API a device is reached through, as it appears in its identifier.
     */
    std::string api_name(sycl::backend backend) {
        switch (backend) {
            case sycl::backend::cuda:       return "cuda";
            case sycl::backend::hip:        return "hip";
            case sycl::backend::level_zero: return "level_zero";
            case sycl::backend::ocl:        return "opencl";
            case sycl::backend::omp:        return host_api;
            default:                        return "other";
        }
    }

    /**
     * @brief Whether @a device is AdaptiveCpp's OpenMP host device, which the kernels are JIT-compiled for through LLVM like any other.
     */
    bool is_host(const sycl::device& device) {
        return device.get_backend() == sycl::backend::omp;
    }

    /**
     * @brief How likely a device reached through @a backend is to be a dedicated GPU; lower is more likely.
     *
     * Nothing in SYCL says whether a GPU is integrated, so the API is used as a stand-in: CUDA and HIP only reach dedicated cards in practice,
     * while integrated GPUs are only reached through Level Zero and OpenCL. Level Zero is preferred over OpenCL, as it is the one Intel's
     * dedicated cards are meant to be used through.
     */
    int rank(sycl::backend backend) {
        switch (backend) {
            case sycl::backend::cuda:
            case sycl::backend::hip:        return 0;
            case sycl::backend::level_zero: return 1;
            case sycl::backend::ocl:        return 2;
            case sycl::backend::omp:        return 4; // the host device, always last
            default:                        return 3;
        }
    }

    /**
     * @brief Every GPU the runtime can reach, followed by the host device, ranked as described for device_count in the ABI header.
     *        Listed once, on first use.
     *
     * All non-ACPP CPU devices are skipped. 
     */
    const std::vector<Candidate>& candidates() {
        static const std::vector<Candidate> list = [] {
            std::vector<Candidate> result;
            try {
                for (const auto& device : sycl::device::get_devices()) {
                    if (!device.is_gpu() && !is_host(device)) {continue;}
                    result.push_back(Candidate{
                        .device=device,
                        .id=api_name(device.get_backend()) + ": " + device.get_info<sycl::info::device::name>()
                    });
                }
            } catch (...) {
                return std::vector<Candidate>{};
            }

            // stable, so that devices of equal rank keep the order the runtime listed them in, and with it their suffixes below
            std::ranges::stable_sort(result, [] (const Candidate& a, const Candidate& b) {
                return rank(a.device.get_backend()) < rank(b.device.get_backend());
            });

            // tell identical devices apart; the suffix is derived from the plain id, so it is checked against those
            std::vector<std::string> plain;
            for (auto& candidate : result) {
                const auto duplicates = std::ranges::count(plain, candidate.id);
                plain.push_back(candidate.id);
                if (duplicates != 0) {candidate.id += " #" + std::to_string(duplicates + 1);}
            }
            return result;
        }();
        return list;
    }

    // index into candidates() of the device to use, as chosen by select_device; the first one if never chosen
    std::size_t chosen_device = 0;

    // whether select_device has been called; the host device is only ever used if it was chosen this way
    bool device_chosen = false;

    // set as soon as the context is first asked for, after which the device can no longer be changed
    bool device_fixed = false;

    // why there is no context, if there is none
    std::string context_error;

    /**
     * @brief The device and the buffers shared by all calculations.
     *
     * Acquiring a device and compiling the kernels is expensive, so both are done once and kept for the lifetime of the process.
     */
    class Context {
        public:
            static Context* get() {
                static std::unique_ptr<Context> instance = [] () -> std::unique_ptr<Context> {
                    device_fixed = true;
                    try {
                        const auto& list = candidates();
                        if (list.empty() || (!device_chosen && is_host(list[chosen_device].device))) {
                            context_error = "no GPU was found";
                            return nullptr;
                        }
                        return std::unique_ptr<Context>(new Context(list[chosen_device].device));
                    } catch (const std::exception& e) {
                        context_error = "could not use the GPU: " + std::string(e.what());
                        return nullptr;
                    } catch (...) {
                        context_error = "could not use the GPU";
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
            std::unordered_map<const float*, Upload> uploaded;              // host -> device, deduplicated
            std::vector<sycl::event> pending;                               // uploads the next kernel must wait for

            /**
             * @brief The host-side sources of the job and tile uploads.
             *
             * A memcpy only promises to have read its source once its event completes; some runtimes (e.g. AdaptiveCpp's OpenMP and OpenCL
             * backends) read it later, from a worker thread. The sources are therefore kept until the next begin(), which waits for the queue.
             * Moving a vector in keeps its buffer, so the pointers handed to memcpy stay valid.
             */
            std::vector<std::vector<DeviceJob>> staged_jobs;
            std::vector<std::vector<TileRef>> staged_tiles;
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
            void* histograms = nullptr; // bin_t<weighted> elements
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

                const std::size_t bytes = static_cast<std::size_t>(slots)*bin_count*bin_bytes(weighted);
                void* fresh = sycl::malloc_device(bytes, queue);
                if (fresh == nullptr) {throw std::runtime_error("sycl_backend: out of device memory");}
                queue.memset(fresh, 0, bytes).wait(); // all-zero bits are 0.0 as well

                if (histograms != nullptr) {
                    queue.wait(); // in-flight kernels are still accumulating into the old buffer
                    const std::size_t old = static_cast<std::size_t>(slot_capacity)*bin_count*bin_bytes(weighted);
                    queue.memcpy(fresh, histograms, old).wait();
                    sycl::free(histograms, queue);
                }
                histograms = fresh;
                slot_capacity = slots;
                zeroing = {}; // everything above was waited for, and the buffer it referred to is gone
            }

        private:
            explicit Context(const sycl::device& device) : queue(device) {
                name = device.get_info<sycl::info::device::name>();
                local_mem = queue.get_device().get_info<sycl::info::device::local_mem_size>();
                warmup();
            }

            // Compile all kernel variants on a trivial input, so the first real calculation does not pay for it.
            void warmup();
    };

    /**
     * @brief Upload any coordinate set this calculation has not seen yet, and return the uploads in job order.
     *
     * The deduplication map lives for the whole calculation rather than one submit, so a coordinate set shared by jobs in different submits 
     * is still uploaded only once - which is the common case for the symmetry managers, where one body's coordinates appear in many jobs. 
     * The uploads are asynchronous; their events are collected so the kernels of this submit can wait on them without the host having to.
     */
    std::vector<std::pair<Upload, Upload>> upload_coordinates(Context& context, const Job* jobs, int n_jobs) {
        auto upload = [&] (const float* host, std::uint32_t n) -> Upload {
            if (auto it = context.uploaded.find(host); it != context.uploaded.end()) {return it->second;}
            const std::size_t bytes = static_cast<std::size_t>(n)*sizeof(sycl::float4);
            auto* device = static_cast<sycl::float4*>(context.arena.allocate(context.queue, bytes));
            context.pending.push_back(context.queue.memcpy(device, host, bytes));

            float max_weight = 0;
            for (std::uint32_t i = 0; i < n; ++i) {max_weight = sycl::max(max_weight, sycl::fabs(host[4*i + 3]));}
            return context.uploaded[host] = Upload{.device=device, .max_weight=max_weight};
        };

        std::vector<std::pair<Upload, Upload>> result(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const Upload a1 = upload(jobs[i].a1, jobs[i].n1);
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
        bin_t<weighted>* histograms, float inv_width, int bin_count, int local_bins, // NOLINT
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

        // only distances in the local bins are summed in fixed point, so they bound the center terms
        const float center_scale = fixed_scale((context.local_bins + 0.5)/context.inv_width);

        const auto coordinates = upload_coordinates(context, jobs, n_jobs);
        std::vector<DeviceJob> device_jobs(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const bool self = jobs[i].a2 == nullptr;
            const auto& [a1, a2] = coordinates[i];
            const float value_scale = fixed_scale(static_cast<double>(a1.max_weight)*a2.max_weight);
            const auto scaling = static_cast<double>(jobs[i].scaling);
            device_jobs[i] = DeviceJob{
                .a1=a1.device, .a2=a2.device, .n1=jobs[i].n1, .n2=self ? jobs[i].n1 : jobs[i].n2,
                .scaling=jobs[i].scaling, .slot=jobs[i].slot, .value_scale=value_scale, .center_scale=center_scale,
                .value_factor=scaling/value_scale, .center_factor=scaling/center_scale
            };
        }

        std::vector<TileRef> tiles, diagonal;
        build_tiles(jobs, n_jobs, tiles, diagonal);
        const std::size_t n_regular = tiles.size();
        tiles.insert(tiles.end(), diagonal.begin(), diagonal.end()); // uploaded as one, regular tiles first

        // every buffer below is arena memory belonging to this submit alone, so the kernels may
        // still be reading it when the next submit arrives
        const std::size_t job_bytes = device_jobs.size()*sizeof(DeviceJob);
        const std::size_t tile_bytes = tiles.size()*sizeof(TileRef);
        auto* device_job_array = static_cast<DeviceJob*>(context.arena.allocate(context.queue, job_bytes));
        auto* device_tiles = static_cast<TileRef*>(context.arena.allocate(context.queue, tile_bytes));

        context.pending.push_back(context.queue.memcpy(device_job_array, device_jobs.data(), job_bytes));
        context.pending.push_back(context.queue.memcpy(device_tiles, tiles.data(), tile_bytes)); // never empty, since n_jobs > 0
        context.staged_jobs.push_back(std::move(device_jobs));
        context.staged_tiles.push_back(std::move(tiles));

        auto* histograms = static_cast<bin_t<weighted>*>(context.histograms);
        launch<weighted, false>(
            context, device_tiles, n_regular, device_job_array,
            histograms, context.inv_width, context.bin_count, context.local_bins, context.pending
        );
        launch<weighted, true>(
            context, device_tiles + n_regular, diagonal.size(), device_job_array,
            histograms, context.inv_width, context.bin_count, context.local_bins, context.pending
        );

        // the uploads are now sequenced before the kernels that read them; a later submit brings its
        // own, and must not make its kernels wait on this one's as well
        context.pending.clear();
    }

    /**
     * @brief Wait for everything submitted since begin(), then copy the histograms straight into @a out, whose layout they already share.
     */
    template<bool weighted>
    void finish_on(Context& context, int n_slots, bin_t<weighted>* out) {
        // check_finish has already bounded n_slots by what was submitted, so nothing has to be allocated or cleared here to answer it
        context.queue.wait_and_throw();

        const std::size_t bins = static_cast<std::size_t>(n_slots)*context.bin_count;
        if (bins != 0) {context.queue.memcpy(out, context.histograms, bins*sizeof(bin_t<weighted>)).wait();}
        context.active = false;
    }

    //  Start a calculation on an idle device, discarding whatever the last one left behind.
    void begin_on(Context& context, int bin_count, float inv_width, bool weighted) {
        context.queue.wait_and_throw(); // nothing may still be reading the arena we are about to rewind
        context.arena.rewind();
        context.uploaded.clear();
        context.pending.clear();
        context.staged_jobs.clear();
        context.staged_tiles.clear();

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
        const std::size_t bytes = static_cast<std::size_t>(context.dirty_slots)*bin_count*bin_bytes(weighted);
        context.zeroing = bytes == 0 ? sycl::event{} : context.queue.memset(context.histograms, 0, bytes);
        context.dirty_slots = 0;
    }

    void Context::warmup() {
        // one self and one cross job, so both the diagonal and the regular kernel are compiled
        const std::vector<float> atoms(4*tile_size, 1.f);
        const std::array<Job, 2> jobs = {
            Job{.a1=atoms.data(), .a2=nullptr, .n1=tile_size, .n2=0, .scaling=2, .slot=0},
            Job{.a1=atoms.data(), .a2=atoms.data(), .n1=tile_size, .n2=tile_size, .scaling=2, .slot=0}
        };
        begin_on(*this, 8, 1.f, false);
        submit_on<false>(*this, jobs.data(), 2);
        std::vector<double> out_unweighted(8);
        finish_on<false>(*this, 1, out_unweighted.data());

        begin_on(*this, 8, 1.f, true);
        submit_on<true>(*this, jobs.data(), 2);
        std::vector<WeightedBin> out_weighted(8);
        finish_on<true>(*this, 1, out_weighted.data());
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
                error_message = "sycl_backend: no usable device (" + context_error + ")";
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

std::int32_t ausaxs::gpu::abi::ausaxs_gpu_device_count() {
    return static_cast<std::int32_t>(candidates().size());
}

const char* ausaxs::gpu::abi::ausaxs_gpu_device_id(std::int32_t index) {
    const auto& list = candidates();
    if (index < 0 || static_cast<std::size_t>(index) >= list.size()) {return nullptr;}
    return list[index].id.c_str();
}

Status ausaxs::gpu::abi::ausaxs_gpu_select_device(const char* id) {
    error_message.clear();
    const auto& list = candidates();
    const auto found = std::ranges::find_if(list, [id] (const Candidate& c) {return id != nullptr && c.id == id;});
    if (found == list.end()) {
        error_message = "sycl_backend: no device is called \"" + std::string(id != nullptr ? id : "") + "\"";
        return Status::invalid_input;
    }

    const auto index = static_cast<std::size_t>(found - list.begin());
    if (device_fixed && index != chosen_device) {
        error_message = "sycl_backend: cannot switch to \"" + found->id + "\", since \"" + list[chosen_device].id + "\" is already in use";
        return Status::invalid_input;
    }
    chosen_device = index;
    device_chosen = true;
    return Status::ok;
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
    return guarded([&] (Context& context) {finish_on<false>(context, n_slots, out);});
}

Status ausaxs::gpu::abi::ausaxs_gpu_finish_weighted(std::int32_t n_slots, WeightedBin* out) {
    if (auto status = check_finish(Context::get(), true, n_slots); status != Status::ok) {return status;}
    return guarded([&] (Context& context) {finish_on<true>(context, n_slots, out);});
}
