// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <gpu/WebGPU/WebGPUBackend.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <utility/observer_ptr.h>

#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief WebGPU kernel for simple histogram calculations. See gpu::WebGPUBackend for the implementation.
     *
     * Submitting a calculation to the GPU costs a few milliseconds of buffer setup and readback latency
     * regardless of its size, which is more than the CPU needs in total for anything but large inputs.
     * The queued calculations are therefore only recorded here, and the backend to run them on is picked
     * in run() from the total number of pair distances they amount to. The caller has to keep the data
     * alive until run() either way, so nothing is lost by deciding that late.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class WebGPUSimple : public SimpleKernel<weighted_bins, variable_bin_width> {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        public:
            using run_result = typename SimpleKernel<weighted_bins, variable_bin_width>::run_result;

            /**
             * @brief The number of pair distances from which the GPU is faster than the CPU.
             *
             * This is a rough break-even point rather than a sharp one; it depends on both the GPU and
             * the number of CPU cores available. Measured on a Radeon RX 9070 XT against 12 CPU threads,
             * the two are within ~20% of each other anywhere between 1e7 and 5e7 pair distances.
             */
            constexpr static std::size_t gpu_threshold = 20'000'000;

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
                        result = replay<Backend>();
                    } catch (const std::exception& e) {
                        // a machine without a usable device should still get its histogram
                        logging::log(std::string("WebGPUSimple::run: falling back to the cpu kernel: ") + e.what());
                        console::print_warning(std::string("WebGPUSimple::run: could not use the GPU: ") + e.what());
                        result = replay<SimpleCPU<weighted_bins, variable_bin_width>>();
                    }
                } else {
                    result = replay<SimpleCPU<weighted_bins, variable_bin_width>>();
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

            // thin adapter giving the backend the kernel interface, so run() can treat both the same way
            struct Backend {
                gpu::WebGPUBackend<weighted_bins, variable_bin_width> backend;
                void enqueue_calculate_self(const CompactCoordinates_t& a, int scaling, int merge_id) {
                    backend.submit_self(a, scaling, merge_id);
                }
                void enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int scaling, int merge_id) {
                    backend.submit_cross(a1, a2, scaling, merge_id);
                }
                run_result run() {return backend.run();}
            };

            template<typename kernel_t>
            run_result replay() {
                kernel_t kernel;
                for (const auto& job : self_jobs) {
                    kernel.enqueue_calculate_self(*job.a1, job.scaling, job.merge_id);
                }
                for (const auto& job : cross_jobs) {
                    kernel.enqueue_calculate_cross(*job.a1, *job.a2, job.scaling, job.merge_id);
                }
                return kernel.run();
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
