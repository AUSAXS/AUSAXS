// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/detail/HistogramManagerMTFFHelpers.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <container/ThreadLocalWrapper.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>
#include <utility/Logging.h>

#include <type_traits>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::container;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool wb, bool vbw>
HistogramManagerMTFFAvg<wb, vbw>::~HistogramManagerMTFFAvg() = default;

template<bool wb, bool vbw>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFAvg<wb, vbw>::calculate() {return calculate_all();}

template<bool wb, bool vbw>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFAvg<wb, vbw>::calculate_all() {
    assert(this->protein != nullptr && "HistogramManagerMTFFAvg::calculate_all: Molecule is not set.");
    logging::log("HistogramManagerMTFFAvg::calculate: starting calculation");

    using GenericDistribution1D_t = typename GenericDistribution1D<wb>::type;
    using GenericDistribution2D_t = typename GenericDistribution2D<wb>::type;
    using GenericDistribution3D_t = typename GenericDistribution3D<wb>::type;
    auto pool = utility::multi_threading::get_global_pool();

    data_a_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_bodies());
    data_w_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_waters());
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<GenericDistribution3D_t> p_aa_all(settings::form_factor::max_ff_types, settings::form_factor::max_ff_types, settings::axes::bin_count); // ff_type1, ff_type2, distance
    auto calc_aa = [&data_a, &p_aa_all, data_a_size] (int imin, int imax) {
        auto& p_aa = p_aa_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = i+1;                    // atom
            for (; j+15 < data_a_size; j+=16) {
                evaluate_aa16<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j+7 < data_a_size; j+=8) {
                evaluate_aa8<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate_aa4<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate_aa1<vbw, 2>(p_aa, data_a, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution2D_t> p_aw_all(settings::form_factor::max_ff_types, settings::axes::bin_count); // ff_type, distance
    auto calc_aw = [&data_w, &data_a, &p_aw_all, data_w_size] (int imin, int imax) {
        auto& p_aw = p_aw_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = 0;                      // water
            for (; j+15 < data_w_size; j+=16) {
                evaluate_aw16<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j+7 < data_w_size; j+=8) {
                evaluate_aw8<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate_aw4<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate_aw1<vbw, 1>(p_aw, data_a, data_w, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_ww_all(settings::axes::bin_count); // distance
    auto calc_ww = [&data_w, &p_ww_all, data_w_size] (int imin, int imax) {
        auto& p_ww = p_ww_all.get();
        for (int i = imin; i < imax; ++i) { // water
            int j = i+1;                    // water
            for (; j+15 < data_w_size; j+=16) {
                evaluate16<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+7 < data_w_size; j+=8) {
                evaluate8<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<vbw, 2>(p_ww, data_w, data_w, i, j);
            }
        }
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    int job_size = settings::general::detail::job_size;
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pool->detach_task(
            [&calc_aa, i, job_size, data_a_size] () {calc_aa(i, std::min(i+job_size, data_a_size));}
        );
    }
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pool->detach_task(
            [&calc_aw, i, job_size, data_a_size] () {calc_aw(i, std::min(i+job_size, data_a_size));}
        );
    }
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        pool->detach_task(
            [&calc_ww, i, job_size, data_w_size] () {calc_ww(i, std::min(i+job_size, data_w_size));}
        );
    }

    pool->wait();
    auto p_aa = p_aa_all.merge();
    auto p_aw = p_aw_all.merge();
    auto p_ww = p_ww_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    for (int i = 0; i < data_a_size; ++i) {
        if constexpr (wb) {
            p_aa.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0, 0.0f);
        } else {
            p_aa.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0);
        }
    }
    if constexpr (wb) {
        p_ww.add_index(0, WeightedEntry(data_w_size, data_w_size, 0));
    } else {
        p_ww.add_index(0, data_w_size);
    }

    GenericDistribution1D_t p_tot(settings::axes::bin_count);
    {   // sum all elements to the total
        unsigned int n_active = form_factor::get_active_count();
        for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < n_active; ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<>());
            }
        }
        for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_aw.begin(ff1), p_tot.begin(), std::plus<>());
        }
        std::transform(p_tot.begin(), p_tot.end(), p_ww.begin(), p_tot.begin(), std::plus<>());
    }

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (unsigned int i = p_tot.size()-1; i >= 10; --i) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    pool->detach_task([&p_aa, max_bin] () { p_aa.resize(max_bin); });
    pool->detach_task([&p_aw, max_bin] () { p_aw.resize(max_bin); });
    pool->detach_task([&p_ww, max_bin] () { p_ww.resize(max_bin); });
    pool->detach_task([&p_tot, max_bin] () { p_tot.resize(max_bin); });
    pool->wait();

    const int n_active = static_cast<int>(form_factor::get_active_count());
    const int ff_start = form_factor::start_index_for_explicit_exv();
    const int nbins = static_cast<int>(max_bin);
    using aa_entry_t = std::remove_cvref_t<decltype(p_aa.index(0, 0, 0))>;

    std::vector<aa_entry_t> xx(nbins, aa_entry_t{});
    for (int a = ff_start; a < n_active; ++a) {
        std::vector<aa_entry_t> ax(nbins, aa_entry_t{});
        for (int b = ff_start; b < n_active; ++b) {
            for (int d = 0; d < nbins; ++d) {
                ax[d] += p_aa.index(a, b, d);
                ax[d] += p_aa.index(b, a, d);
                xx[d] += p_aa.index(a, b, d);
            }
        }
        for (int d = 0; d < nbins; ++d) {
            p_aa.index(a, form_factor::exv_bin, d) = 0.5*ax[d];
        }
    }
    for (int d = 0; d < nbins; ++d) {
        p_aa.index(form_factor::exv_bin, form_factor::exv_bin, d) = xx[d];
    }

    for (int a = ff_start; a < n_active; ++a) {
        for (int d = 0; d < nbins; ++d) {
            p_aw.index(form_factor::exv_bin, d) += p_aw.index(a, d);
        }
    }

    // multiply the excluded volume charge onto the excluded volume bins
    double Z_exv_avg = 
        this->protein->size_atom() == 0 
        ? 0 
        : this->protein->get_volume_grid()*constants::charge::density::water/this->protein->size_atom()
    ;
    for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < form_factor::get_active_count(); ++ff1) {
        std::transform(p_aa.begin(ff1, form_factor::exv_bin), p_aa.end(ff1, form_factor::exv_bin), p_aa.begin(ff1, form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg;});
    }
    std::transform(p_aa.begin(form_factor::exv_bin, form_factor::exv_bin), p_aa.end(form_factor::exv_bin, form_factor::exv_bin), p_aa.begin(form_factor::exv_bin, form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg*Z_exv_avg;});
    std::transform(p_aw.begin(form_factor::exv_bin), p_aw.end(form_factor::exv_bin), p_aw.begin(form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg;});

    return std::make_unique<CompositeDistanceHistogramFFAvg>(
        std::move(Distribution3D(std::move(p_aa))), 
        std::move(Distribution2D(std::move(p_aw))), 
        std::move(Distribution1D(std::move(p_ww))), 
        std::move(p_tot)
    );
}

template class hist::HistogramManagerMTFFAvg<false, false>;
template class hist::HistogramManagerMTFFAvg<false, true>;
template class hist::HistogramManagerMTFFAvg<true, false>;
template class hist::HistogramManagerMTFFAvg<true, true>;