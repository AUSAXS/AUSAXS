#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <utility/MultiThreading.h>
#include <container/ThreadLocalWrapper.h>

#include <cassert>

using namespace ausaxs;

// This class could be improved by dispatching the calculations as soon as they are enqueued. 
// The primary difficulty with implementing this is to ensure there is enough space in the results vector, 
// as dynamically extending it invalidates the references stored in the lambdas. 

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int merge_id) {
    self.emplace_back(a);
    int res_id = merge_id == -1 ? self.size()-1 : merge_id;
    self_merge_ids.emplace_back(res_id);
    return res_id;
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int merge_id) {
    cross_1.emplace_back(a1);
    cross_2.emplace_back(a2);
    int res_id = merge_id == -1 ? cross_1.size()-1 : merge_id;
    cross_merge_ids.emplace_back(res_id);
    return res_id;
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::size_self_result() const {
    return self.size();
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::size_cross_result() const {
    assert(cross_1.size() == cross_2.size() && "Cross-correlation data size mismatch");
    return cross_1.size();
}

template<bool weighted_bins>
hist::distance_calculator::SimpleCalculator<weighted_bins>::run_result hist::distance_calculator::SimpleCalculator<weighted_bins>::run() {
    auto pool = utility::multi_threading::get_global_pool();
    assert(cross_1.size() == cross_2.size() && "Cross-correlation data size mismatch");

    // calculate self-correlations
    std::vector<container::ThreadLocalWrapper<GenericDistribution1D_t>> results_self(
        self_merge_ids.back()+1, 
        container::ThreadLocalWrapper<GenericDistribution1D_t>(constants::axes::d_axis.bins)
    );
    for (unsigned int data_idx = 0; data_idx < self.size(); ++data_idx) {
        int res_idx = self_merge_ids[data_idx];
        int data_a_size = (int) self[res_idx].get().size();
        int job_size = settings::general::detail::job_size;

        // calculate upper triangle
        for (int i = 0; i < data_a_size; i+=job_size) {
            pool->detach_task(
                [this, &results_self, res_idx, data_idx, imin = i, imax = std::min(i+job_size, data_a_size)] () {
                    int data_a_size = static_cast<int>(self[data_idx].get().size());
                    auto& p_aa = results_self[res_idx].get();
                    auto& data_a = self[data_idx].get();
                    for (int i = imin; i < imax; ++i) { // atom
                        int j = i+1;                    // atom
                        for (; j+7 < data_a_size; j+=8) {
                            evaluate8<weighted_bins, 2>(p_aa, data_a, data_a, i, j);
                        }

                        for (; j+3 < data_a_size; j+=4) {
                            evaluate4<weighted_bins, 2>(p_aa, data_a, data_a, i, j);
                        }

                        for (; j < data_a_size; ++j) {
                            evaluate1<weighted_bins, 2>(p_aa, data_a, data_a, i, j);
                        }
                    }
                }
            );
        }

        // calculate skipped diagonal
        pool->detach_task(
            [this, &results_self, res_idx, data_idx] () {
                auto& p_aa = results_self[res_idx].get();
                auto& data_a = self[data_idx].get();
                p_aa.add(0, std::accumulate(
                    data_a.get_data().begin(), 
                    data_a.get_data().end(), 
                    0.0, 
                    [] (double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}
                ));
            }
        );
    }

    // calculate cross-correlations
    std::vector<container::ThreadLocalWrapper<GenericDistribution1D_t>> results_cross(
        cross_merge_ids.back()+1, 
        container::ThreadLocalWrapper<GenericDistribution1D_t>(constants::axes::d_axis.bins)
    );
    for (unsigned int data_idx = 0; data_idx < cross_1.size(); ++data_idx) {
        int res_idx = cross_merge_ids[data_idx];
        int data_b_size = static_cast<int>(cross_2[data_idx].get().size());
        int job_size = settings::general::detail::job_size;
        
        for (int i = 0; i < data_b_size; i+=job_size) {
            pool->detach_task(
                [this, &results_cross, res_idx, data_idx, imin = i, imax = std::min(i+job_size, data_b_size)] () {
                    int data_a_size = static_cast<int>(cross_1[data_idx].get().size());
                    auto& p_ab = results_cross[res_idx].get();
                    auto& data_a = cross_1[data_idx].get();
                    auto& data_b = cross_2[data_idx].get();
                    for (int i = imin; i < imax; ++i) { // b
                        int j = 0;                      // a
                        for (; j+7 < data_a_size; j+=8) {
                            evaluate8<weighted_bins, 2>(p_ab, data_b, data_a, i, j);
                        }

                        for (; j+3 < data_a_size; j+=4) {
                            evaluate4<weighted_bins, 2>(p_ab, data_b, data_a, i, j);
                        }

                        for (; j < data_a_size; ++j) {
                            evaluate1<weighted_bins, 2>(p_ab, data_b, data_a, i, j);
                        }
                    }
                }
            );
        }
    }

    pool->wait();
    run_result result(self_merge_ids.back(), cross_merge_ids.back());

    // results_X contains the thread-local histograms. We need to merge them into a single final histogram
    for (int i = 0; i < static_cast<int>(self.size()); ++i) {
        result.self[self_merge_ids[i]] = results_self[i].merge();
    }

    for (int i = 0; i < static_cast<int>(cross_1.size()); ++i) {
        result.cross[cross_merge_ids[i]] = results_cross[i].merge();
    }

    // cleanup
    self.clear();
    cross_1.clear();
    cross_2.clear();

    return result;
}

template class hist::distance_calculator::SimpleCalculator<true>;
template class hist::distance_calculator::SimpleCalculator<false>;