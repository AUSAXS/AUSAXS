// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/GPULoader.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distance_calculator/detail/CalculatorCPU.h>
#include <utility/observer_ptr.h>

#include <cassert>
#include <cstdint>
#include <deque>
#include <memory>
#include <span>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator::detail {
    /**
     * @brief Simple histogram calculation on whichever GPU backend is installed, into the rows of a HistogramStore.
     *        If the device fails, or if there is none, the CPU calculator is used instead.
     *
     * @tparam unit_weights Whether every point weighs 1. The device always multiplies the weights, so they are sent as 1.
     */
    template<bool weighted_bins, bool variable_bin_width, bool unit_weights>
    class GPUKernel {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using Row = std::span<typename HistogramStore<weighted_bins>::entry_type>;
        public:
            /**
             * @brief Construct a kernel accumulating into @a store, which must outlive it.
             */
            explicit GPUKernel(HistogramStore<weighted_bins>& store) : store(&store) {
                if (!gpu::GPULoader::available()) {switch_to_cpu();}
            }

            void enqueue_calculate_self(const CompactCoordinates_t& a, Row row, int scaling) {
                if (a.empty()) {store->reset(row); return;}
                open_session();
                if (on_cpu) {cpu->enqueue_calculate_self(a, row, scaling); return;}

                jobs.emplace_back(Job{&a, nullptr, row, scaling});
                int slot = resolve(row);
                diagonal[slot] += self_weight(a, scaling);
                submit(gpu::abi::Job{
                    coordinates(a), nullptr, static_cast<std::uint32_t>(a.size()), 0,
                    static_cast<std::uint32_t>(2*scaling), static_cast<std::uint32_t>(slot)
                });
            }

            void enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, Row row, int pair_factor) {
                if (a1.empty() || a2.empty()) {store->reset(row); return;}
                open_session();
                if (on_cpu) {cpu->enqueue_calculate_cross(a1, a2, row, pair_factor); return;}

                jobs.emplace_back(Job{&a1, &a2, row, pair_factor});
                int slot = resolve(row);
                submit(gpu::abi::Job{
                    coordinates(a1), coordinates(a2),
                    static_cast<std::uint32_t>(a1.size()), static_cast<std::uint32_t>(a2.size()),
                    static_cast<std::uint32_t>(pair_factor), static_cast<std::uint32_t>(slot)
                });
            }

            void hold() {
                // holds do not nest: an inner release_hold() would dispatch the outer group early
                assert(!holding && "GPUKernel::hold: already holding");
                holding = true;
            }

            void release_hold() {
                assert(holding && "GPUKernel::release_hold: not holding");
                holding = false;
                flush();
            }

            void run() {
                // sanity check: the caller should always remember to release a held group
                assert(queued.empty() && "GPUKernel::run: the held group was never released");
                flush();

                // the cpu calculator folds by itself; the device writes its rows directly, but the resets still need to be folded
                if (on_cpu) {cpu->run();}
                else {
                    store->fold();
                    read_back();
                }

                // cleanup
                jobs.clear();
                slots.clear();
                rows.clear();
                diagonal.clear();
                coordinate_buffers.clear();
                queued.clear();
                session_open = false;
                holding = false;
            }

        private:
            struct Job {
                observer_ptr<const CompactCoordinates_t> a1, a2; // a2 is null for self-correlations
                Row row;
                int factor; // the scaling of a self-correlation, or the pair factor of a cross-correlation
            };

            observer_ptr<HistogramStore<weighted_bins>> store;
            std::vector<gpu::abi::Job> queued;                      // held jobs, dispatched by release_hold()
            std::vector<Job> jobs;                                  // kept only to replay on the cpu if the device fails
            std::unordered_map<const void*, int> slots;             // row -> the device histogram it accumulates into, keyed by its first bin
            std::vector<Row> rows;                                  // device histogram -> the row it belongs to
            std::vector<double> diagonal;                           // per device histogram, the zero-distance contribution
            bool session_open = false;                              // whether begin() has been issued for the batch being built
            bool holding = false;                                   // whether jobs are being collected into a group, see hold()
            std::deque<std::vector<float>> coordinate_buffers;
            std::unique_ptr<CalculatorCPU<weighted_bins, variable_bin_width, unit_weights>> cpu;
            bool on_cpu = false;                                    // whether the device was given up on, see switch_to_cpu()

            /**
             * @brief Open a device session, on the first job of a batch. Does nothing on later jobs.
             *        A device that refuses to start hands this batch, and every later one, to the CPU calculator.
             */
            void open_session() {
                if (session_open || on_cpu) {return;}
                session_open = true;

                const auto& backend = gpu::GPULoader::get();
                auto status = backend.begin(
                    static_cast<std::int32_t>(store->bins()),
                    hist::detail::WidthController<variable_bin_width>::get_inv_width(),
                    weighted_bins
                );
                check_status(status, "open a device session");
            }

            /**
             * @brief Give up on the device and replay whatever was already submitted to it.
             *
             * Anything still queued on the device is simply abandoned; the next begin() waits for it before reusing the memory it holds. 
             */
            void switch_to_cpu() {
                cpu = std::make_unique<CalculatorCPU<weighted_bins, variable_bin_width, unit_weights>>(*store);
                on_cpu = true;
                for (const auto& job : jobs) {
                    if (job.a2 == nullptr) {cpu->enqueue_calculate_self(*job.a1, job.row, job.factor);}
                    else {cpu->enqueue_calculate_cross(*job.a1, *job.a2, job.row, job.factor);}
                }
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
             * @brief Hand the calculation to the CPU calculator unless the device call succeeded.
             * @param action What was attempted, for the warning GPULoader prints. See report_failure().
             */
            void check_status(gpu::abi::Status status, std::string_view action) {
                if (status == gpu::abi::Status::ok) {return;}
                gpu::GPULoader::report_failure(action, status);
                switch_to_cpu();
            }

            /**
             * @brief Wait for the device and assign each returned histogram, plus its diagonal, to its row.
             */
            void read_back() {
                if (!session_open) {return;}

                const int bin_count = store->bins();
                const int n_slots = static_cast<int>(rows.size());
                const auto& backend = gpu::GPULoader::get();

                if constexpr (weighted_bins) {
                    std::vector<gpu::abi::WeightedBin> out(static_cast<std::size_t>(n_slots)*bin_count);
                    auto status = backend.finish_weighted(n_slots, out.data());
                    if (status != gpu::abi::Status::ok) {replay_on_cpu(status); return;}
                    for (int slot = 0; slot < n_slots; ++slot) {
                        auto row = rows[slot];
                        for (int i = 0; i < bin_count; ++i) {
                            const auto& bin = out[static_cast<std::size_t>(slot)*bin_count + i];
                            row[i] = hist::detail::WeightedEntry{bin.value, bin.count, bin.center};
                        }
                        row[0] += hist::detail::WeightedEntry{diagonal[slot], static_cast<std::int64_t>(diagonal[slot]), 0};
                    }
                } else {
                    std::vector<double> out(static_cast<std::size_t>(n_slots)*bin_count);
                    auto status = backend.finish_unweighted(n_slots, out.data());
                    if (status != gpu::abi::Status::ok) {replay_on_cpu(status); return;}
                    for (int slot = 0; slot < n_slots; ++slot) {
                        auto row = rows[slot];
                        std::copy_n(out.begin() + static_cast<std::ptrdiff_t>(slot)*bin_count, bin_count, row.begin());
                        row[0] += diagonal[slot];
                    }
                }
            }

            /**
             * @brief A device that failed at the last moment should still produce a histogram.
             */
            void replay_on_cpu(gpu::abi::Status status) {
                gpu::GPULoader::report_failure("read the histograms back", status);
                switch_to_cpu();
                cpu->run();
            }

            const float* coordinates(const CompactCoordinates_t& a) {
                coordinate_buffers.emplace_back(a.size()*4);
                auto& packed = coordinate_buffers.back();
                for (int i = 0; i < a.size(); ++i) {
                    packed[4*i] = a.x(i);
                    packed[4*i + 1] = a.y(i);
                    packed[4*i + 2] = a.z(i);
                    packed[4*i + 3] = unit_weights ? 1 : a.get_weight(i);
                }
                return packed.data();
            }

            /**
             * @brief The contribution of the zero distance of every atom with itself.
             */
            static double self_weight(const CompactCoordinates_t& a, int scaling) {
                if constexpr (unit_weights) {return static_cast<double>(scaling)*a.size();}
                double total_weight = 0;
                for (int i = 0; i < a.size(); ++i) {
                    double weight = a.get_weight(i);
                    total_weight += weight*weight;
                }
                return scaling*total_weight;
            }

            /**
             * @brief The device histogram of @a row, allocating one the first time it is seen in this batch.
             *        Only rows with device work get one, so every histogram finish() is asked for was submitted to.
             */
            int resolve(Row row) {
                if (auto it = slots.find(row.data()); it != slots.end()) {return it->second;}
                int slot = static_cast<int>(rows.size());
                rows.push_back(row);
                diagonal.push_back(0);
                return slots[row.data()] = slot;
            }
    };
}
