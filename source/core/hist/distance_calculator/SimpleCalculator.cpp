#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <utility/MultiThreading.h>

#include <cassert>

using namespace ausaxs;

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& data, int merge_id) {
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
                        evaluate8<weighted_bins, 2>(p_aa, data, data, i, j);
                    }

                    for (; j+3 < data_size; j+=4) {
                        evaluate4<weighted_bins, 2>(p_aa, data, data, i, j);
                    }

                    for (; j < data_size; ++j) {
                        evaluate1<weighted_bins, 2>(p_aa, data, data, i, j);
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

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_cross(
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
                        evaluate8<weighted_bins, 2>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j+3 < data_1_size; j+=4) {
                        evaluate4<weighted_bins, 2>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j < data_1_size; ++j) {
                        evaluate1<weighted_bins, 2>(p_ab, data_2, data_1, i, j);
                    }
                }
            }
        );
    }

    return res_idx;
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::size_self_result() const {
    return self_results.size();
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::size_cross_result() const {
    return cross_results.size();
}

template<bool weighted_bins>
hist::distance_calculator::SimpleCalculator<weighted_bins>::run_result hist::distance_calculator::SimpleCalculator<weighted_bins>::run() {
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

template class hist::distance_calculator::SimpleCalculator<true>;
template class hist::distance_calculator::SimpleCalculator<false>;