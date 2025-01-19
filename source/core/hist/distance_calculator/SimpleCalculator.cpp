#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <utility/MultiThreading.h>
#include <container/ThreadLocalWrapper.h>

#include <cassert>

using namespace ausaxs;

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& a) {
    self.push_back(a);
    return self.size()-1;
}

template<bool weighted_bins>
int hist::distance_calculator::SimpleCalculator<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2) {
    cross_1.push_back(a1);
    cross_2.push_back(a2);
    return cross_1.size()-1;
}

template<bool weighted_bins>
hist::distance_calculator::SimpleCalculator<weighted_bins>::run_result hist::distance_calculator::SimpleCalculator<weighted_bins>::run() {
    auto pool = utility::multi_threading::get_global_pool();
    assert(cross_1.size() == cross_2.size() && "Cross-correlation data size mismatch");

    // calculate self-correlations
    std::vector<container::ThreadLocalWrapper<GenericDistribution1D_t>> results_self(self.size(), container::ThreadLocalWrapper<GenericDistribution1D_t>(constants::axes::d_axis.bins));
    for (unsigned int idx = 0; idx < self.size(); ++idx) {
        int data_a_size = (int) self[idx].get().size();
        int job_size = settings::general::detail::job_size;

        // calculate upper triangle
        for (int i = 0; i < data_a_size; i+=job_size) {
            pool->detach_task(
                [this, &results_self, idx, imin = i, imax = std::min(i+job_size, data_a_size)] () {
                    int data_a_size = (int) self[idx].get().size();
                    auto& p_aa = results_self[idx].get();
                    auto& data_a = self[idx].get();
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
            [this, &results_self, idx] () {
                auto& p_aa = results_self[idx].get();
                auto& data_a = self[idx].get();
                p_aa.add(0, std::accumulate(
                    data_a.get_data().begin(), 
                    data_a.get_data().end(), 
                    0.0, 
                    [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}
                ));
            }
        );
    }

    // calculate cross-correlations
    std::vector<container::ThreadLocalWrapper<GenericDistribution1D_t>> results_cross(cross_1.size(), container::ThreadLocalWrapper<GenericDistribution1D_t>(constants::axes::d_axis.bins));
    std::vector<std::function<void(int, int)>> tasks_cross(cross_1.size()); // to avoid unnecessary copies
    for (unsigned int idx = 0; idx < cross_1.size(); ++idx) {
        int data_b_size = (int) cross_2[idx].get().size();
        int job_size = settings::general::detail::job_size;
        
        for (int i = 0; i < data_b_size; i+=job_size) {
            pool->detach_task(
                [this, &results_cross, idx, imin = i, imax = std::min(i+job_size, data_b_size)] () {
                    int data_a_size = (int) cross_1[idx].get().size();
                    auto& p_ab = results_cross[idx].get();
                    auto& data_a = cross_1[idx].get();
                    auto& data_b = cross_2[idx].get();
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
    run_result result(self.size(), cross_1.size());

    // results_X contains the thread-local histograms. We need to merge them into a single final histogram
    for (int i = 0; i < static_cast<int>(self.size()); ++i) {
        result.self[i] = results_self[i].merge();
    }

    for (int i = 0; i < static_cast<int>(cross_1.size()); ++i) {
        result.cross[i] = results_cross[i].merge();
    }

    return result;
}

template class hist::distance_calculator::SimpleCalculator<true>;
template class hist::distance_calculator::SimpleCalculator<false>;