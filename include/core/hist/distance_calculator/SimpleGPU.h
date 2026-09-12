// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/GPULoader.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/observer_ptr.h>

#include <cassert>
#include <cstdint>
#include <memory>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Simple histogram calculation on whichever GPU backend is installed.
     *        If the device fails, or if there is none, the CPU kernel is used instead.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SimpleGPU {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        public:
            using run_result = typename SimpleCPU<weighted_bins, variable_bin_width>::run_result;

            /**
             * @brief Construct a kernel whose result histograms span @a bin_count bins.
             */
            explicit SimpleGPU(int bin_count) : bin_count(bin_count) {
                if (!gpu::GPULoader::available()) {switch_to_cpu();}
            }

            int enqueue_calculate_self(const CompactCoordinates_t& a, int scaling = 1, int merge_id = -1) {
                open_session();
                if (on_cpu) {return cpu->enqueue_calculate_self(a, scaling, merge_id);}

                auto [slot, index] = resolve(self_slots, merge_id);
                self_jobs.emplace_back(Job{&a, nullptr, scaling, merge_id});
                diagonal[slot] += self_weight(a, scaling);
                submit(gpu::abi::Job{
                    coordinates(a), nullptr, static_cast<std::uint32_t>(a.size()), 0,
                    static_cast<std::uint32_t>(scaling), static_cast<std::uint32_t>(slot)
                });
                return index;
            }

            int enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int scaling = 1, int merge_id = -1) {
                open_session();
                if (on_cpu) {return cpu->enqueue_calculate_cross(a1, a2, scaling, merge_id);}

                auto [slot, index] = resolve(cross_slots, merge_id);
                cross_jobs.emplace_back(Job{&a1, &a2, scaling, merge_id});
                submit(gpu::abi::Job{
                    coordinates(a1), coordinates(a2),
                    static_cast<std::uint32_t>(a1.size()), static_cast<std::uint32_t>(a2.size()),
                    static_cast<std::uint32_t>(scaling), static_cast<std::uint32_t>(slot)
                });
                return index;
            }

            int size_self_result() const {
                if (on_cpu) {return cpu->size_self_result();}
                return static_cast<int>(self_slots.size());
            }

            int size_cross_result() const {
                if (on_cpu) {return cpu->size_cross_result();}
                return static_cast<int>(cross_slots.size());
            }

            void hold() {holding = true;}

            void release_hold() {
                holding = false;
                flush();
            }

            run_result run() {
                // sanity check: the caller should always remember to release a held group
                assert(queued.empty() && "SimpleGPU::run: the held group was never released");
                flush();

                run_result result = on_cpu ? cpu->run() : read_back();

                // cleanup
                self_jobs.clear();
                cross_jobs.clear();
                self_slots.clear();
                cross_slots.clear();
                diagonal.clear();
                queued.clear();
                next_slot = 0;
                session_open = false;
                holding = false;

                return result;
            }

        private:
            struct Job {
                observer_ptr<const CompactCoordinates_t> a1, a2; // a2 is null for self-correlations
                int scaling;
                int merge_id; // as resolved by resolve(), never -1
            };

            struct Slot {
                int slot;   // index of the device histogram, from one counter shared by self and cross
                int index;  // position among the results of its own kind, which is what callers index by
            };

            int bin_count;                                          // bins spanned by every histogram in this batch
            std::vector<gpu::abi::Job> queued;                      // held jobs, dispatched by release_hold()
            std::vector<Job> self_jobs{}, cross_jobs{};                 // kept only to replay on the cpu if the device fails
            std::unordered_map<int, Slot> self_slots{}, cross_slots{};  // merge id -> where its result is
            bool session_open = false;                              // whether begin() has been issued for the batch being built
            bool holding = false;                                   // whether jobs are being collected into a group, see hold()
            std::vector<double> diagonal;                           // per slot, the zero-distance contribution
            int next_slot = 0;
            std::unique_ptr<SimpleCPU<weighted_bins, variable_bin_width>> cpu;
            bool on_cpu = false;                                    // whether the device was given up on, see switch_to_cpu()

            /**
             * @brief Open a device session, on the first job of a batch. Does nothing on later jobs.
             *        A device that refuses to start hands this batch, and every later one, to the cpu kernel.
             */
            void open_session() {
                if (session_open || on_cpu) {return;}
                session_open = true;

                const auto& backend = gpu::GPULoader::get();
                auto status = backend.begin(
                    static_cast<std::int32_t>(bin_count),
                    hist::detail::WidthController<variable_bin_width>::get_inv_width(),
                    weighted_bins
                );
                check_status(status, "open a device session");
            }

            /**
             * @brief Give up on the device and replay whatever was already submitted to it.
             *
             * Anything still queued on the device is simply abandoned; the next begin() waits for it
             * before reusing the memory it holds.
             */
            void switch_to_cpu() {
                cpu = std::make_unique<SimpleCPU<weighted_bins, variable_bin_width>>(bin_count);
                on_cpu = true;
                for (const auto& job : self_jobs) {cpu->enqueue_calculate_self(*job.a1, job.scaling, job.merge_id);}
                for (const auto& job : cross_jobs) {cpu->enqueue_calculate_cross(*job.a1, *job.a2, job.scaling, job.merge_id);}
            }

            void submit(const gpu::abi::Job& job) {
                queued.push_back(job);
                if (!holding) {flush();}
            }

            /**
             * @brief Hand everything queued so far to the device as a single submission.
             */
            void flush() {
                if (queued.empty() || on_cpu) {queued.clear(); return;}

                const auto& backend = gpu::GPULoader::get();
                auto status = backend.submit(queued.data(), static_cast<std::int32_t>(queued.size()));
                queued.clear();
                check_status(status, "submit a group of jobs");
            }

            /**
             * @brief Hand the calculation to the cpu kernel unless the device call succeeded.
             * @param action What was attempted, for the warning GPULoader prints. See report_failure().
             */
            void check_status(gpu::abi::Status status, std::string_view action) {
                if (status == gpu::abi::Status::ok) {return;}
                gpu::GPULoader::report_failure(action, status);
                switch_to_cpu();
            }

            /**
             * @brief Wait for the device and turn the returned histograms into distributions.
             *        A batch with no jobs in it never opened a session, and finish() without a matching
             *        begin() is an error, so there is nothing to wait for and nothing to read.
             */
            run_result read_back() {
                if (!session_open) {return run_result{};}

                const int bin_count = this->bin_count;
                const auto& backend = gpu::GPULoader::get();
                std::vector<GenericDistribution1D_t> slots(next_slot);

                if constexpr (weighted_bins) {
                    std::vector<gpu::abi::WeightedBin> out(static_cast<std::size_t>(next_slot)*bin_count);
                    auto status = backend.finish_weighted(next_slot, out.data());
                    if (status != gpu::abi::Status::ok) {return replay_on_cpu(status);}
                    for (int slot = 0; slot < next_slot; ++slot) {
                        slots[slot] = GenericDistribution1D_t(bin_count);
                        for (int i = 0; i < bin_count; ++i) {
                            const auto& bin = out[static_cast<std::size_t>(slot)*bin_count + i];
                            slots[slot].add_index(i, hist::detail::WeightedEntry{
                                bin.value,
                                bin.count,
                                bin.center
                            });
                        }
                        if (diagonal[slot] == 0) {continue;}
                        slots[slot].add_index(0, hist::detail::WeightedEntry{
                            diagonal[slot], static_cast<std::int64_t>(diagonal[slot]), 0
                        });
                    }
                } else {
                    std::vector<double> out(static_cast<std::size_t>(next_slot)*bin_count);
                    auto status = backend.finish_unweighted(next_slot, out.data());
                    if (status != gpu::abi::Status::ok) {return replay_on_cpu(status);}
                    for (int slot = 0; slot < next_slot; ++slot) {
                        slots[slot] = GenericDistribution1D_t(bin_count);
                        for (int i = 0; i < bin_count; ++i) {
                            slots[slot].add_index(i, out[static_cast<std::size_t>(slot)*bin_count + i]);
                        }
                        if (diagonal[slot] != 0) {slots[slot].add_index(0, diagonal[slot]);}
                    }
                }

                run_result result;
                for (const auto& [merge_id, where] : self_slots) {result.self[merge_id] = slots[where.slot];}
                for (const auto& [merge_id, where] : cross_slots) {result.cross[merge_id] = slots[where.slot];}
                return result;
            }

            /**
             * @brief A device that failed at the last moment should still produce a histogram.
             */ 
            run_result replay_on_cpu(gpu::abi::Status status) {
                gpu::GPULoader::report_failure("read the histograms back", status);
                switch_to_cpu();
                return cpu->run();
            }

            static const float* coordinates(const CompactCoordinates_t& a) {
                static_assert(
                    sizeof(typename std::decay_t<decltype(a.get_data())>::value_type) == 4*sizeof(float),
                    "The coordinates must be tightly packed [x, y, z, w] floats to match the kernel layout."
                );
                return reinterpret_cast<const float*>(a.get_data().data());
            }

            /**
             * @brief The contribution of the zero distance of every atom with itself.
             */
            static double self_weight(const CompactCoordinates_t& a, int scaling) {
                return scaling*std::accumulate(
                    a.get_data().begin(), a.get_data().end(), 0.0,
                    [] (double sum, const auto& val) {return sum + val.value.w*val.value.w;}
                );
            }

            /**
             * @brief Assign a merge id its slot, allocating one the first time it is seen.
             */
            Slot resolve(std::unordered_map<int, Slot>& slots, int& merge_id) {
                if (auto it = slots.find(merge_id); merge_id != -1 && it != slots.end()) {return it->second;}
                if (merge_id == -1) {merge_id = static_cast<int>(slots.size());}

                Slot where{next_slot++, static_cast<int>(slots.size())};
                diagonal.resize(next_slot, 0);
                return slots[merge_id] = where;
            }
    };
}
