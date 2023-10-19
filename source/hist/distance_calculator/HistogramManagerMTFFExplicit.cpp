#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/CompositeDistanceHistogramFFExplicit.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactor.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <container/Container3D.h>
#include <container/Container2D.h>
#include <container/Container1D.h>
#include <constants/Constants.h>

#include <BS_thread_pool.hpp>

using namespace container;
using namespace hist;

HistogramManagerMTFFExplicit::HistogramManagerMTFFExplicit(HistogramManager& hm) : HistogramManager(hm) {}

HistogramManagerMTFFExplicit::~HistogramManagerMTFFExplicit() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit::calculate() {return calculate_all();}

std::unique_ptr<CompositeDistanceHistogram> HistogramManagerMTFFExplicit::calculate_all() {
    hist::detail::CompactCoordinatesFF data_p(protein->get_bodies());
    hist::detail::CompactCoordinatesFF data_h = hist::detail::CompactCoordinatesFF(protein->get_waters());

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    BS::thread_pool pool(settings::general::threads);
    auto calc_pp = [&data_p] (unsigned int imin, unsigned int imax) {
        Container3D<double> p_aa(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
        Container3D<double> p_ax(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
        Container3D<double> p_xx(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = i+1;
            for (; j+7 < data_p.get_size(); j+=8) {
                auto res = data_p[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {
                    p_aa.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2*res.weight[k];
                    p_ax.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2*data_p[i+k].value.w;
                    p_xx.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2;
                }
            }

            for (; j+3 < data_p.get_size(); j+=4) {
                auto res = data_p[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {
                    p_aa.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2*res.weight[k];
                    p_ax.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2*data_p[i+k].value.w;
                    p_xx.index(data_p.get_ff_type(i+k), data_p.get_ff_type(j+k), res.distance[k]) += 2;
                }
            }

            for (; j < data_p.get_size(); ++j) {
                auto res = data_p[i].evaluate(data_p[j]);
                p_aa.index(data_p.get_ff_type(i), data_p.get_ff_type(j), res.distance) += 2*res.weight;
                p_ax.index(data_p.get_ff_type(i), data_p.get_ff_type(j), res.distance) += 2*data_p[i].value.w;
                p_xx.index(data_p.get_ff_type(i), data_p.get_ff_type(j), res.distance) += 2;
             }
        }
        return std::tuple(std::move(p_aa), std::move(p_ax), std::move(p_xx));
    };

    auto calc_hp = [&data_h, &data_p] (unsigned int imin, unsigned int imax) {
        Container2D<double> p_ww(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
        Container2D<double> p_wx(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = 0;
            for (; j+7 < data_p.get_size(); j+=8) {
                auto res = data_h[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3], data_p[j+4], data_p[j+5], data_p[j+6], data_p[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {
                    p_ww.index(data_p.get_ff_type(j+k), res.distance[k]) += res.weight[k];
                    p_wx.index(data_p.get_ff_type(j+k), res.distance[k]) += 2;
                }
            }

            for (; j+3 < data_p.get_size(); j+=4) {
                auto res = data_h[i].evaluate(data_p[j], data_p[j+1], data_p[j+2], data_p[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {
                    p_ww.index(data_p.get_ff_type(j+k), res.distance[k]) += res.weight[k];
                    p_wx.index(data_p.get_ff_type(j+k), res.distance[k]) += 2;
                }
            }

            for (; j < data_p.get_size(); ++j) {
                auto res = data_h[i].evaluate(data_p[j]);
                p_ww.index(data_p.get_ff_type(j), res.distance) += res.weight;
                p_wx.index(data_p.get_ff_type(j), res.distance) += 2;
            }
        }
        return std::tuple(std::move(p_ww), std::move(p_wx));
    };

    auto calc_hh = [&data_h] (unsigned int imin, unsigned int imax) {
        Container1D<double> p_hh(constants::axes::d_axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = i+1;
            for (; j+7 < data_h.get_size(); j+=8) {
                auto res = data_h[i].evaluate(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3], data_h[j+4], data_h[j+5], data_h[j+6], data_h[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {
                    p_hh.index(res.distance[k]) += 2*res.weight[k];
                }
            }

            for (; j+3 < data_h.get_size(); j+=4) {
                auto res = data_h[i].evaluate(data_h[j], data_h[j+1], data_h[j+2], data_h[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {
                    p_hh.index(res.distance[k]) += 2*res.weight[k];
                }
            }

            for (; j < data_h.get_size(); ++j) {
                auto res = data_h[i].evaluate(data_h[j]);
                p_hh.index(res.distance) += 2*res.weight;
            }
        }
        return p_hh;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    unsigned int job_size = settings::general::detail::job_size;
    BS::multi_future<std::tuple<Container3D<double>, Container3D<double>, Container3D<double>>> pp;
    for (unsigned int i = 0; i < data_p.get_size(); i+=job_size) {
        pp.push_back(pool.submit(calc_pp, i, std::min(i+job_size, data_p.get_size())));
    }
    BS::multi_future<std::tuple<Container2D<double>, Container2D<double>>> hp;
    for (unsigned int i = 0; i < data_h.get_size(); i+=job_size) {
        hp.push_back(pool.submit(calc_hp, i, std::min(i+job_size, data_h.get_size())));
    }
    BS::multi_future<Container1D<double>> hh;
    for (unsigned int i = 0; i < data_h.get_size(); i+=job_size) {
        hh.push_back(pool.submit(calc_hh, i, std::min(i+job_size, data_h.get_size())));
    }

    //#################//
    // COLLECT RESULTS //
    //#################//
    auto p_pp_future = pool.submit(
        [&]() {
            Container3D<double> p_aa(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            Container3D<double> p_ax(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            Container3D<double> p_xx(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            for (const auto& [p_aa_tmp, p_ax_tmp, p_xx_tmp] : pp.get()) {
                std::transform(p_aa.begin(), p_aa.end(), p_aa_tmp.begin(), p_aa.begin(), std::plus<double>());
                std::transform(p_ax.begin(), p_ax.end(), p_ax_tmp.begin(), p_ax.begin(), std::plus<double>());
                std::transform(p_xx.begin(), p_xx.end(), p_xx_tmp.begin(), p_xx.begin(), std::plus<double>());
            }
            return std::tuple(std::move(p_aa), std::move(p_ax), std::move(p_xx));
        }
    );

    auto p_hp_future = pool.submit(
        [&]() {
            Container2D<double> p_wa(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
            Container2D<double> p_wx(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
            for (const auto& [p_wa_tmp, p_wx_tmp] : hp.get()) {
                std::transform(p_wa.begin(), p_wa.end(), p_wa_tmp.begin(), p_wa.begin(), std::plus<double>());
                std::transform(p_wx.begin(), p_wx.end(), p_wx_tmp.begin(), p_wx.begin(), std::plus<double>());
            }
            return std::tuple(std::move(p_wa), std::move(p_wx));
        }
    );

    auto p_hh_future = pool.submit(
        [&]() {
            Container1D<double> p_ww(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hh.get()) {
                std::transform(p_ww.begin(), p_ww.end(), tmp.begin(), p_ww.begin(), std::plus<double>());
            }
            return p_ww;
        }
    );
    pool.wait_for_tasks();
    auto [p_aa, p_ax, p_xx] = p_pp_future.get();
    auto [p_wa, p_wx] = p_hp_future.get();
    auto p_ww = p_hh_future.get();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    for (unsigned int i = 0; i < data_p.get_size(); ++i) {p_aa.index(data_p.get_ff_type(i), data_p.get_ff_type(i), 0) += std::pow(data_p[i].value.w, 2);}
    p_ww.index(0) = std::accumulate(data_h.get_data().begin(), data_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& data) {return sum + std::pow(data.value.w, 2);});

    // this is counter-intuitive, but splitting the loop into separate parts is likely faster since it allows both SIMD optimizations and better cache usage
    std::vector<double> p_tot(constants::axes::d_axis.bins, 0);
    {   // sum all elements to the total
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<double>());
            }
        }
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_wa.begin(ff1), p_tot.begin(), std::plus<double>());
        }
        std::transform(p_tot.begin(), p_tot.end(), p_ww.begin(), p_tot.begin(), std::plus<double>());
    }

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (unsigned int i = p_tot.size()-1; i >= 10; --i) {
        if (p_tot[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    Container3D<double> p_aa_short(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), max_bin);
    Container3D<double> p_ax_short(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), max_bin);
    Container3D<double> p_xx_short(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), max_bin);
    Container2D<double> p_wa_short(form_factor::get_count_without_excluded_volume(), max_bin);
    Container2D<double> p_wx_short(form_factor::get_count_without_excluded_volume(), max_bin);
    Container1D<double> p_ww_short(max_bin);
    {   // move first max_bin elements to a new container
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                std::move(p_aa.begin(ff1, ff2), p_aa.begin(ff1, ff2)+max_bin, p_aa_short.begin(ff1, ff2));
                std::move(p_ax.begin(ff1, ff2), p_ax.begin(ff1, ff2)+max_bin, p_ax_short.begin(ff1, ff2));
                std::move(p_xx.begin(ff1, ff2), p_xx.begin(ff1, ff2)+max_bin, p_xx_short.begin(ff1, ff2));
            }
            std::move(p_wa.begin(ff1), p_wa.begin(ff1)+max_bin, p_wa_short.begin(ff1));
            std::move(p_wx.begin(ff1), p_wx.begin(ff1)+max_bin, p_wx_short.begin(ff1));
        }
        std::move(p_ww.begin(), p_ww.begin()+max_bin, p_ww_short.begin());
    }
    p_tot.resize(max_bin);
    return std::make_unique<CompositeDistanceHistogramFFExplicit>(std::move(p_aa_short), std::move(p_ax_short), std::move(p_xx_short),
                                                                  std::move(p_wa_short), std::move(p_wx_short), std::move(p_ww_short),
                                                                  std::move(p_tot), Axis(0, max_bin*constants::axes::d_axis.width(), max_bin));
}