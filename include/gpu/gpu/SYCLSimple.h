// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <hist/detail/data/WidthControllers.h>
#include <gpu/SYCL/SYCLBackend.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <utility/observer_ptr.h>

#include <numeric>
#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief SYCL kernel for simple histogram calculations. See gpu::sycl_backend for the implementation.
     *
     * Like the WebGPU kernel, submitting a calculation to the GPU costs a fixed amount of setup and
     * readback latency that is more than the CPU needs in total for small inputs, so the queued
     * calculations are only recorded here and the backend to run them on is picked in run().
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SYCLSimple : public SimpleKernel<weighted_bins, variable_bin_width> {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        public:
            using run_result = typename SimpleKernel<weighted_bins, variable_bin_width>::run_result;

            /**
             * @brief The number of pair distances from which the GPU is faster than the CPU.
             *
             * The SYCL kernel dispatches through the native driver rather than a portability layer, so
             * its fixed cost is a good deal lower than the WebGPU one and it pays off much earlier.
             */
            constexpr static std::size_t gpu_threshold = 2'000'000;

            int enqueue_calculate_self(const CompactCoordinates_t& a, int scaling = 1, int merge_id = -1) override {
                int index = resolve(self_merge_ids, merge_id);
                self_jobs.emplace_back(Job{&a, nullptr, scaling, merge_id});
                pairs += a.size()*(a.size() - 1)/2;
                return index;
            }

            int enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int scaling = 1, int merge_id = -1) override {
                int index = resolve(cross_merge_ids, merge_id);
                cross_jobs.emplace_back(Job{&a1, &a2, scaling, merge_id});
                pairs += a1.size()*a2.size();
                return index;
            }

            int size_self_result() const override {return static_cast<int>(self_merge_ids.size());}
            int size_cross_result() const override {return static_cast<int>(cross_merge_ids.size());}

            run_result run() override {
                run_result result;
                if (gpu_threshold <= pairs) {
                    try {
                        result = run_gpu();
                    } catch (const std::exception& e) {
                        // a machine without a usable device should still get its histogram
                        logging::log(std::string("SYCLSimple::run: falling back to the cpu kernel: ") + e.what());
                        console::print_warning(std::string("SYCLSimple::run: could not use the GPU: ") + e.what());
                        result = run_cpu();
                    }
                } else {
                    result = run_cpu();
                }
                self_jobs.clear();
                cross_jobs.clear();
                self_merge_ids.clear();
                cross_merge_ids.clear();
                pairs = 0;
                return result;
            }

        private:
            struct Job {
                observer_ptr<const CompactCoordinates_t> a1, a2; // a2 is null for self-correlations
                int scaling;
                int merge_id; // as resolved by resolve(), never -1
            };

            std::vector<Job> self_jobs, cross_jobs;
            std::unordered_map<int, int> self_merge_ids, cross_merge_ids;
            std::size_t pairs = 0;

            static const float* coordinates(const CompactCoordinates_t& a) {
                static_assert(
                    sizeof(typename std::decay_t<decltype(a.get_data())>::value_type) == 4*sizeof(float),
                    "The coordinates must be tightly packed [x, y, z, w] floats to match the kernel layout."
                );
                return reinterpret_cast<const float*>(a.get_data().data());
            }

            /**
             * @brief The contribution of the zero distance of every atom with itself.
             *
             * The kernel only evaluates the pairs above the diagonal, so this is added here exactly as
             * the CPU kernel does it — as a single term per job, rather than one per atom, so the two
             * agree to the last bit.
             */
            static double self_weight(const CompactCoordinates_t& a, int scaling) {
                return scaling*std::accumulate(
                    a.get_data().begin(), a.get_data().end(), 0.0,
                    [] (double sum, const auto& val) {return sum + val.value.w*val.value.w;}
                );
            }

            run_result run_cpu() {
                SimpleCPU<weighted_bins, variable_bin_width> kernel;
                for (const auto& job : self_jobs) {
                    kernel.enqueue_calculate_self(*job.a1, job.scaling, job.merge_id);
                }
                for (const auto& job : cross_jobs) {
                    kernel.enqueue_calculate_cross(*job.a1, *job.a2, job.scaling, job.merge_id);
                }
                return kernel.run();
            }

            run_result run_gpu() {
                const int bin_count = static_cast<int>(settings::axes::bin_count);
                const int n_self = static_cast<int>(self_merge_ids.size());
                const int n_slots = n_self + static_cast<int>(cross_merge_ids.size());

                // the self results take the first slots, the cross results the rest
                std::vector<gpu::sycl_backend::Job> jobs;
                jobs.reserve(self_jobs.size() + cross_jobs.size());
                std::vector<double> diagonal(n_slots, 0);
                for (const auto& job : self_jobs) {
                    const int slot = self_merge_ids.at(job.merge_id);
                    jobs.emplace_back(gpu::sycl_backend::Job{
                        coordinates(*job.a1), nullptr, static_cast<std::uint32_t>(job.a1->size()), 0,
                        static_cast<std::uint32_t>(job.scaling), static_cast<std::uint32_t>(slot)
                    });
                    diagonal[slot] += self_weight(*job.a1, job.scaling);
                }
                for (const auto& job : cross_jobs) {
                    jobs.emplace_back(gpu::sycl_backend::Job{
                        coordinates(*job.a1), coordinates(*job.a2),
                        static_cast<std::uint32_t>(job.a1->size()), static_cast<std::uint32_t>(job.a2->size()),
                        static_cast<std::uint32_t>(job.scaling),
                        static_cast<std::uint32_t>(n_self + cross_merge_ids.at(job.merge_id))
                    });
                }

                const float inv_width = hist::detail::WidthController<variable_bin_width>::get_inv_width();
                std::vector<GenericDistribution1D_t> slots(n_slots);
                if constexpr (weighted_bins) {
                    std::vector<gpu::sycl_backend::WeightedBin> out(static_cast<std::size_t>(n_slots)*bin_count);
                    gpu::sycl_backend::run_weighted(jobs.data(), static_cast<int>(jobs.size()), inv_width, bin_count, n_slots, out.data());
                    for (int slot = 0; slot < n_slots; ++slot) {
                        slots[slot] = GenericDistribution1D_t(bin_count);
                        for (int i = 0; i < bin_count; ++i) {
                            const auto& bin = out[static_cast<std::size_t>(slot)*bin_count + i];
                            slots[slot].add_index(i, hist::detail::WeightedEntry{
                                bin.value/gpu::sycl_backend::fixed_point_scale,
                                static_cast<unsigned int>(bin.count),
                                bin.center/gpu::sycl_backend::fixed_point_scale
                            });
                        }
                    }
                    for (int slot = 0; slot < n_slots; ++slot) {
                        if (diagonal[slot] == 0) {continue;}
                        slots[slot].add_index(0, hist::detail::WeightedEntry{
                            diagonal[slot], static_cast<unsigned int>(diagonal[slot]), 0
                        });
                    }
                } else {
                    std::vector<std::int64_t> out(static_cast<std::size_t>(n_slots)*bin_count);
                    gpu::sycl_backend::run_unweighted(jobs.data(), static_cast<int>(jobs.size()), inv_width, bin_count, n_slots, out.data());
                    for (int slot = 0; slot < n_slots; ++slot) {
                        slots[slot] = GenericDistribution1D_t(bin_count);
                        for (int i = 0; i < bin_count; ++i) {
                            slots[slot].add_index(i, out[static_cast<std::size_t>(slot)*bin_count + i]/gpu::sycl_backend::fixed_point_scale);
                        }
                        if (diagonal[slot] != 0) {slots[slot].add_index(0, diagonal[slot]);}
                    }
                }

                run_result result;
                for (const auto& [merge_id, index] : self_merge_ids) {
                    result.self[merge_id] = slots[index];
                }
                for (const auto& [merge_id, index] : cross_merge_ids) {
                    result.cross[merge_id] = slots[n_self + index];
                }
                return result;
            }

            /**
             * @brief Assign a job to a result slot exactly like the kernels themselves would, so the
             *        returned indices and the result keys are the same whichever backend ends up running.
             *        The merge_id of the job is replaced by the resolved one.
             */
            int resolve(std::unordered_map<int, int>& merge_ids, int& merge_id) {
                if (auto it = merge_ids.find(merge_id); merge_id != -1 && it != merge_ids.end()) {
                    return it->second;
                }
                int index = static_cast<int>(merge_ids.size());
                if (merge_id == -1) {merge_id = index;}
                merge_ids[merge_id] = index;
                return index;
            }
    };
}
