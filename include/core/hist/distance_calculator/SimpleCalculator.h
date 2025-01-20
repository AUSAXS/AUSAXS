#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <container/ThreadLocalWrapper.h>
#include <utility/MultiThreading.h>

#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Simple interface to queue histogram calculations. 
     *        Submit the data and then call calculate to get the result.
     *        The caller must guarantee the lifetime of all submitted data.
     */
    template<bool weighted_bins>
    class SimpleCalculator {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        struct run_result {
            run_result(int size_self, int size_cross) : self(size_self), cross(size_cross) {}
            std::vector<GenericDistribution1D_t> self;
            std::vector<GenericDistribution1D_t> cross;
        };

        public:
            /**
             * @brief Queue a self-correlation calculation. 
             *        This is faster than calling the cross-correlation method with the same data, as some optimizations can be made. 
             *
             * @param a The data to calculate the self-correlation for. The reference must be valid until calculate is called.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources. 
             * @param scaling The scaling factor to apply to the result.
             *
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling = 1, int merge_id = -1);

            /**
             * @brief Queue a cross-correlation calculation. 
             *
             * @param a1 The first set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param a2 The second set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources.
             * @param scaling The scaling factor to apply to the result.
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling = 1, int merge_id = -1);

            /**
             * @brief Get the current size of the result vector. 
             */
            int size_self_result() const; 
            int size_cross_result() const; //< @copydoc size_self_result

            /**
             * @brief Calculate the queued histograms. 
             *        This will block until all calculations are done.
             *
             * @return The calculated histograms. 
             */
            run_result run();

        private:
            std::vector<std::unique_ptr<container::ThreadLocalWrapper<GenericDistribution1D_t>>> self_results, cross_results;

            template<int scaling>
            int enqueue_calculate_self(const hist::detail::CompactCoordinates& data, int merge_id);

            template<int scaling>
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates& data_1, const hist::detail::CompactCoordinates& data_2, int merge_id);
    };
}

template<bool weighted_bins> template<int scaling>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_self(
    const hist::detail::CompactCoordinates& data, 
    int merge_id
) {
    auto pool = utility::multi_threading::get_global_pool();

    int res_idx;
    if (merge_id == -1 || merge_id == static_cast<int>(self_results.size())) {
        // second condition is to enable using size_self_result() to get the next merge_id
        self_results.emplace_back(std::make_unique<container::ThreadLocalWrapper<GenericDistribution1D_t>>(constants::axes::d_axis.bins));
        res_idx = self_results.size()-1;
    } else {
        assert(merge_id < static_cast<int>(self_results.size()));
        res_idx = merge_id;
    }

    auto res_ptr = self_results[res_idx].get();
    int data_size = static_cast<int>(data.size());
    int job_size = settings::general::detail::job_size;

    // calculate upper triangle
    for (int i = 0; i < data_size; i+=job_size) {
        pool->detach_task(
            [&data, res_ptr, data_size, imin = i, imax = std::min(i+job_size, data_size)] () {
                auto& p_aa = res_ptr->get();
                for (int i = imin; i < imax; ++i) { // atom
                    int j = i+1;                    // atom
                    for (; j+7 < data_size; j+=8) {
                        evaluate8<weighted_bins, 2*scaling>(p_aa, data, data, i, j);
                    }

                    for (; j+3 < data_size; j+=4) {
                        evaluate4<weighted_bins, 2*scaling>(p_aa, data, data, i, j);
                    }

                    for (; j < data_size; ++j) {
                        evaluate1<weighted_bins, 2*scaling>(p_aa, data, data, i, j);
                    }
                }
            }
        );
    }

    // calculate skipped diagonal
    pool->detach_task(
        [&data, res_ptr] () {
            auto& p_aa = res_ptr->get();
            p_aa.add(0, std::accumulate(
                data.get_data().begin(), 
                data.get_data().end(), 
                0.0, 
                [] (double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}
            ));
        }
    );

    return res_idx;
}

template<bool weighted_bins> template<int scaling>
int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_cross(
    const hist::detail::CompactCoordinates& data_1, 
    const hist::detail::CompactCoordinates& data_2, 
    int merge_id
) {
    auto pool = utility::multi_threading::get_global_pool();
    int res_idx;
    if (merge_id == -1) {
        cross_results.emplace_back(std::make_unique<container::ThreadLocalWrapper<GenericDistribution1D_t>>(constants::axes::d_axis.bins));
        res_idx = cross_results.size()-1;
    } else {
        res_idx = merge_id;
    }

    auto res_ptr = cross_results[res_idx].get();
    int data_1_size = static_cast<int>(data_1.size());
    int data_2_size = static_cast<int>(data_2.size());
    int job_size = settings::general::detail::job_size;

    for (int i = 0; i < data_2_size; i+=job_size) {
        pool->detach_task(
            [&data_1, &data_2, res_ptr, data_1_size, imin = i, imax = std::min(i+job_size, data_2_size)] () {
                auto& p_ab = res_ptr->get();
                for (int i = imin; i < imax; ++i) { // b
                    int j = 0;                      // a
                    for (; j+7 < data_1_size; j+=8) {
                        evaluate8<weighted_bins, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j+3 < data_1_size; j+=4) {
                        evaluate4<weighted_bins, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j < data_1_size; ++j) {
                        evaluate1<weighted_bins, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }
                }
            }
        );
    }

    return res_idx;
}

template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::size_self_result() const {
    return self_results.size();
}

template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::size_cross_result() const {
    return cross_results.size();
}

template<bool weighted_bins>
inline ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::run_result ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::run() {
    auto pool = utility::multi_threading::get_global_pool();
    pool->wait();
    run_result result(size_self_result(), size_cross_result());

    // results_X contains the thread-local histograms. We need to merge them into a single final histogram
    for (int i = 0; i < static_cast<int>(result.self.size()); ++i) {
        result.self[i] = self_results[i]->merge();
    }

    for (int i = 0; i < static_cast<int>(result.cross.size()); ++i) {
        result.cross[i] = cross_results[i]->merge();
    }

    // cleanup
    self_results.clear();
    cross_results.clear();

    return result;
}

// this approach is not sustainable
// should higher scaling factors be needed, new add1, add4, and add8 functions should be created which accepts the scaling factor as a parameter
// for now, this is primarily meant for rigidbody optimizations, where larger symmetries are not expected
template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_self(
    const hist::detail::CompactCoordinates& data, 
    int scaling,
    int merge_id
) {
    switch (scaling) {
        case 1: return enqueue_calculate_self<1>(data, merge_id);
        case 2: return enqueue_calculate_self<2>(data, merge_id);
        case 3: return enqueue_calculate_self<3>(data, merge_id);
        case 4: return enqueue_calculate_self<4>(data, merge_id);
        case 5: return enqueue_calculate_self<5>(data, merge_id);
        case 6: return enqueue_calculate_self<6>(data, merge_id);
        case 7: return enqueue_calculate_self<7>(data, merge_id);
        case 8: return enqueue_calculate_self<8>(data, merge_id);
        case 9: return enqueue_calculate_self<9>(data, merge_id);
        case 10: return enqueue_calculate_self<10>(data, merge_id);
        default: throw std::runtime_error("SimpleCalculator::enqueue_calculate_self: too large scaling factor");
    }
}

// this approach is not sustainable
// should higher scaling factors be needed, new add1, add4, and add8 functions should be created which accepts the scaling factor as a parameter
// for now, this is primarily meant for rigidbody optimizations, where larger symmetries are not expected
template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_cross(
    const hist::detail::CompactCoordinates& data_1, 
    const hist::detail::CompactCoordinates& data_2, 
    int scaling,
    int merge_id
) {
    switch (scaling) {
        case 1: return enqueue_calculate_cross<1>(data_1, data_2, merge_id);
        case 2: return enqueue_calculate_cross<2>(data_1, data_2, merge_id);
        case 3: return enqueue_calculate_cross<3>(data_1, data_2, merge_id);
        case 4: return enqueue_calculate_cross<4>(data_1, data_2, merge_id);
        case 5: return enqueue_calculate_cross<5>(data_1, data_2, merge_id);
        case 6: return enqueue_calculate_cross<6>(data_1, data_2, merge_id);
        case 7: return enqueue_calculate_cross<7>(data_1, data_2, merge_id);
        case 8: return enqueue_calculate_cross<8>(data_1, data_2, merge_id);
        case 9: return enqueue_calculate_cross<9>(data_1, data_2, merge_id);
        case 10: return enqueue_calculate_cross<10>(data_1, data_2, merge_id);
        default: throw std::runtime_error("SimpleCalculator::enqueue_calculate_cross: too large scaling factor");
    }
}