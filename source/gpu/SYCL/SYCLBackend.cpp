// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// SYCL implementation of the simple distance histogram kernel, mirroring
// hist::distance_calculator::SimpleCPU. This file is compiled by the SYCL compiler rather than the
// compiler used for the rest of the project, so it only communicates through the trivially copyable
// types of SYCLBackend.h.
//
// The kernel is structured around three decisions:
//
// - Bins are 64-bit fixed-point integers. Accumulating in f32 is far too imprecise, as the bins
//   reach values around 1e11 where the single-precision spacing is several thousand. Unlike WebGPU,
//   which has to emulate this in two 32-bit words, the device atomics here are natively 64 bit.
//
// - Contributions go into a workgroup-local histogram that is flushed once per tile. Only the first
//   local_bins bins fit in workgroup memory; the tail falls back to global memory, which is rare in
//   practice since it covers the first local_bins*bin_width angstrom.
//
// - Work is split into a uniform tile grid over both index ranges, and each thread keeps
//   rows_per_thread atoms of the first range in registers while streaming the second. The streamed
//   atom is the same for every thread in the workgroup, so it is read once into scalar registers and
//   its cost is amortized over all the atoms held in registers. This is what makes the loop
//   arithmetic-bound rather than load-bound.

#include <gpu/SYCL/SYCLBackend.h>

#include <sycl/sycl.hpp>

#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace ausaxs::gpu::sycl_backend;

namespace {
    using i64 = std::int64_t;

    constexpr int workgroup_size = 256;
    constexpr int rows_per_thread = 4;                          // atoms of the first range held in registers
    constexpr int tile_size = workgroup_size*rows_per_thread;   // atoms covered by one tile in either direction

    // as many bins as fit in workgroup memory. the weighted bins need three accumulators each, so
    // fewer of them fit; both leave enough of the 64 kB free for several workgroups per compute unit.
    constexpr int local_bins_unweighted = 2048;
    constexpr int local_bins_weighted = 800;

    constexpr float fixed_scale = static_cast<float>(fixed_point_scale);

    /// @brief Convert a value to its fixed-point representation.
    inline i64 to_fixed(float value) {
        return static_cast<i64>(sycl::rint(value*fixed_scale));
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

    /// @brief A queued correlation as the kernel sees it.
    struct DeviceJob {
        const sycl::float4* a1;
        const sycl::float4* a2; // equal to a1 for self-correlations
        std::uint32_t n1, n2;
        std::uint32_t scaling;
        std::uint32_t slot;
    };

    /// @brief One tile of one job, i.e. the work of a single workgroup.
    struct TileRef {
        std::uint32_t job;
        std::uint32_t i_block;
        std::uint32_t j_block;
    };

    /**
     * @brief Evaluate one tile of the distance matrix into the workgroup-local histogram.
     *
     * @tparam weighted Whether to track the bin count and the summed distance alongside the value.
     * @tparam diagonal Whether the tile lies on the diagonal of a self-correlation, in which case only
     *                  the pairs above the diagonal are evaluated. This is a template parameter so the
     *                  predicate disappears from the inner loop of the (far more common) other tiles.
     */
    template<bool weighted, bool diagonal>
    struct Kernel {
        constexpr static int local_bins = weighted ? local_bins_weighted : local_bins_unweighted;

        const TileRef* tiles;
        const DeviceJob* jobs;
        i64* histograms;        // n_slots*bin_count values, or the same number of {value, count, center} triples
        float inv_width;
        std::uint32_t bin_count;

        // workgroup-local accumulators; count and center are unused when not weighted
        using local_array = sycl::local_accessor<i64, 1>;
        local_array local_value, local_count, local_center;

        void add_to_bin(std::uint32_t bin, float value, float center, std::uint32_t count, i64* global) const {
            if (bin_count <= bin) {return;} // distances beyond the axis are discarded rather than folded into the last bin
            if (bin < static_cast<std::uint32_t>(local_bins)) {
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
            sycl::float4 held[rows_per_thread];
            int index[rows_per_thread];
            bool active[rows_per_thread];
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
                        static_cast<std::uint32_t>(static_cast<int>(sycl::rint(inv_width*d))),
                        w*held[k].w(), scale*d, count, global
                    );
                }
            }

            sycl::group_barrier(group);
            const int bins = sycl::min(local_bins, static_cast<int>(bin_count));
            for (int bin = lid; bin < bins; bin += workgroup_size) {
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

    /// @brief A device allocation that is reused across runs and only ever grows.
    class Scratch {
        public:
            void* get(sycl::queue& queue, std::size_t bytes) {
                if (bytes <= size) {return pointer;}
                release(queue);
                pointer = sycl::malloc_device(bytes, queue);
                if (!pointer) {throw std::runtime_error("sycl_backend: out of device memory");}
                size = bytes;
                return pointer;
            }

            void release(sycl::queue& queue) {
                if (pointer) {sycl::free(pointer, queue);}
                pointer = nullptr;
                size = 0;
            }

        private:
            void* pointer = nullptr;
            std::size_t size = 0;
    };

    /**
     * @brief The device and the buffers shared by all calculations.
     *
     * Acquiring a device and compiling the kernels is expensive, and a rigidbody optimization runs
     * thousands of histograms, so both are done once and kept for the lifetime of the process.
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
            Scratch coordinates, jobs, tiles, output;

        private:
            Context() : queue(sycl::gpu_selector_v) {
                name = queue.get_device().get_info<sycl::info::device::name>();
                warmup();
            }

            /// @brief Compile all kernel variants on a trivial input, so the first real calculation does not pay for it.
            void warmup();
    };

    /// @brief Upload each distinct coordinate set once and return the device pointers, in job order.
    std::vector<std::pair<const sycl::float4*, const sycl::float4*>> upload_coordinates(
        Context& context, const Job* jobs, int n_jobs
    ) {
        std::unordered_map<const float*, std::size_t> offsets; // host pointer -> offset into the device buffer
        std::vector<std::pair<const float*, std::uint32_t>> uploads;
        std::size_t total = 0;
        auto reserve = [&] (const float* host, std::uint32_t n) {
            if (!host || offsets.count(host)) {return;}
            offsets[host] = total;
            uploads.emplace_back(host, n);
            total += n;
        };
        for (int i = 0; i < n_jobs; ++i) {
            reserve(jobs[i].a1, jobs[i].n1);
            if (jobs[i].a2) {reserve(jobs[i].a2, jobs[i].n2);}
        }

        auto* device = static_cast<sycl::float4*>(
            context.coordinates.get(context.queue, sycl::max(total, std::size_t{1})*sizeof(sycl::float4))
        );
        for (const auto& [host, n] : uploads) {
            context.queue.memcpy(device + offsets.at(host), host, static_cast<std::size_t>(n)*sizeof(sycl::float4));
        }

        std::vector<std::pair<const sycl::float4*, const sycl::float4*>> result(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const sycl::float4* a1 = device + offsets.at(jobs[i].a1);
            result[i] = {a1, jobs[i].a2 ? device + offsets.at(jobs[i].a2) : a1};
        }
        return result;
    }

    /**
     * @brief Split the queued jobs into the tiles that make up one workgroup of work each.
     *
     * Self-correlations only enumerate the tiles on and above the diagonal, which are returned
     * separately as they need the extra predicate.
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
                    diagonal.emplace_back(TileRef{job, p, p});
                    for (std::uint32_t q = p + 1; q < n; ++q) {
                        regular.emplace_back(TileRef{job, p, q});
                    }
                }
            } else {
                for (std::uint32_t p = 0; p < blocks(jobs[i].n1); ++p) {
                    for (std::uint32_t q = 0; q < blocks(jobs[i].n2); ++q) {
                        regular.emplace_back(TileRef{job, p, q});
                    }
                }
            }
        }
    }

    template<bool weighted, bool diagonal>
    void launch(
        Context& context, const TileRef* tiles, std::size_t n_tiles, const DeviceJob* jobs,
        i64* histograms, float inv_width, int bin_count
    ) {
        if (n_tiles == 0) {return;}
        constexpr int local_bins = Kernel<weighted, diagonal>::local_bins;
        context.queue.submit([&] (sycl::handler& handler) {
            sycl::local_accessor<i64, 1> value{sycl::range<1>(local_bins), handler};
            sycl::local_accessor<i64, 1> count{sycl::range<1>(weighted ? local_bins : 1), handler};
            sycl::local_accessor<i64, 1> center{sycl::range<1>(weighted ? local_bins : 1), handler};
            handler.parallel_for(
                sycl::nd_range<1>(n_tiles*workgroup_size, workgroup_size),
                Kernel<weighted, diagonal>{
                    tiles, jobs, histograms, inv_width, static_cast<std::uint32_t>(bin_count),
                    value, count, center
                }
            );
        });
    }

    /**
     * @brief Run all queued jobs and write the histograms into @a out.
     *
     * The context is passed in rather than looked up, since warmup() has to be able to call this while
     * the shared instance is still being constructed.
     */
    template<bool weighted>
    void run_on(Context* context, const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, i64* out) {
        if (!context) {throw std::runtime_error("sycl_backend: no usable device");}
        if (n_jobs == 0 || n_slots == 0 || bin_count == 0) {return;}

        constexpr int stride = weighted ? 3 : 1;
        const std::size_t values = static_cast<std::size_t>(n_slots)*bin_count*stride;
        auto* histograms = static_cast<i64*>(context->output.get(context->queue, values*sizeof(i64)));
        context->queue.fill(histograms, i64{0}, values);

        const auto coordinates = upload_coordinates(*context, jobs, n_jobs);
        std::vector<DeviceJob> device_jobs(n_jobs);
        for (int i = 0; i < n_jobs; ++i) {
            const bool self = jobs[i].a2 == nullptr;
            device_jobs[i] = DeviceJob{
                coordinates[i].first, coordinates[i].second, jobs[i].n1,
                self ? jobs[i].n1 : jobs[i].n2, jobs[i].scaling, jobs[i].slot
            };
        }

        std::vector<TileRef> regular, diagonal;
        build_tiles(jobs, n_jobs, regular, diagonal);

        auto* device_job_array = static_cast<DeviceJob*>(
            context->jobs.get(context->queue, device_jobs.size()*sizeof(DeviceJob))
        );
        auto* device_tiles = static_cast<TileRef*>(
            context->tiles.get(context->queue, sycl::max(regular.size() + diagonal.size(), std::size_t{1})*sizeof(TileRef))
        );
        context->queue.memcpy(device_job_array, device_jobs.data(), device_jobs.size()*sizeof(DeviceJob));
        if (!regular.empty()) {
            context->queue.memcpy(device_tiles, regular.data(), regular.size()*sizeof(TileRef));
        }
        if (!diagonal.empty()) {
            context->queue.memcpy(device_tiles + regular.size(), diagonal.data(), diagonal.size()*sizeof(TileRef));
        }
        context->queue.wait(); // the kernels read the buffers filled above

        launch<weighted, false>(*context, device_tiles, regular.size(), device_job_array, histograms, inv_width, bin_count);
        launch<weighted, true>(*context, device_tiles + regular.size(), diagonal.size(), device_job_array, histograms, inv_width, bin_count);
        context->queue.wait_and_throw();
        context->queue.memcpy(out, histograms, values*sizeof(i64)).wait();
    }

    void Context::warmup() {
        // one self and one cross job, so both the diagonal and the regular kernel are compiled
        const std::vector<float> atoms(4*tile_size, 1.f);
        const Job jobs[2] = {
            Job{atoms.data(), nullptr, tile_size, 0, 1, 0},
            Job{atoms.data(), atoms.data(), tile_size, tile_size, 1, 0}
        };
        std::vector<i64> out(3*8);
        run_on<false>(this, jobs, 2, 1.f, 8, 1, out.data());
        run_on<true>(this, jobs, 2, 1.f, 8, 1, out.data());
    }

    template<bool weighted>
    void run(const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, i64* out) {
        run_on<weighted>(Context::get(), jobs, n_jobs, inv_width, bin_count, n_slots, out);
    }
}

bool ausaxs::gpu::sycl_backend::available() {
    return Context::get() != nullptr;
}

const char* ausaxs::gpu::sycl_backend::device_name() {
    Context* context = Context::get();
    return context ? context->name.c_str() : "none";
}

void ausaxs::gpu::sycl_backend::run_unweighted(
    const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, std::int64_t* out
) {
    ::run<false>(jobs, n_jobs, inv_width, bin_count, n_slots, out);
}

void ausaxs::gpu::sycl_backend::run_weighted(
    const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, WeightedBin* out
) {
    static_assert(sizeof(WeightedBin) == 3*sizeof(i64), "WeightedBin must be tightly packed to match the kernel layout.");
    ::run<true>(jobs, n_jobs, inv_width, bin_count, n_slots, reinterpret_cast<i64*>(out));
}
