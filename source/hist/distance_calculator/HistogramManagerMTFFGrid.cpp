#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <hydrate/Grid.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <container/Container2D.h>
#include <container/Container1D.h>
#include <form_factor/FormFactorType.h>

#include <BS_thread_pool.hpp>

using namespace hist;

HistogramManagerMTFFGrid::~HistogramManagerMTFFGrid() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGrid::calculate() {

}

std::unique_ptr<CompositeDistanceHistogram> HistogramManagerMTFFGrid::calculate_all() {
    constexpr unsigned int exv_bin = static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
    hist::detail::CompactCoordinates data_x = hist::detail::CompactCoordinates(protein->get_grid()->generate_excluded_volume());
    auto& data_p = *data_p_ptr;
    auto& data_h = *data_h_ptr;

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    auto calc_xx = [&data_x] (unsigned int imin, unsigned int imax) {
        // container::Container2D<double> p_xx(form_factor::get_count(), constants::axes::d_axis.bins, 0);
        container::Container1D<double> p_xx(constants::axes::d_axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) { // exv
            unsigned int j = 0;                      // exv
            for (; j+7 < data_x.get_size(); j+=8) {
                auto res = data_x[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3], data_x[j+4], data_x[j+5], data_x[j+6], data_x[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_xx.index(res.distance[k]) += 2*res.weight[k];}
            }

            for (; j+3 < data_x.get_size(); j+=4) {
                auto res = data_x[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_xx.index(res.distance[k]) += 2*res.weight[k];}
            }

            for (; j < data_x.get_size(); ++j) {
                auto res = data_x[i].evaluate(data_x[j]);
                p_xx.index(res.distance) += 2*res.weight;
            }
        }
        return p_xx;
    };

    auto calc_ax = [&data_p, &data_x] (unsigned int imin, unsigned int imax) {
        container::Container1D<double> p_ax(constants::axes::d_axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) { // atoms
            unsigned int j = 0;                      // exv
            for (; j+7 < data_x.get_size(); j+=8) {
                auto res = data_p[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3], data_x[j+4], data_x[j+5], data_x[j+6], data_x[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_ax.index(res.distance[k]) += 2*res.weight[k];}
            }

            for (; j+3 < data_x.get_size(); j+=4) {
                auto res = data_p[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_ax.index(res.distance[k]) += 2*res.weight[k];}
            }

            for (; j < data_x.get_size(); ++j) {
                auto res = data_p[i].evaluate(data_x[j]);
                p_ax.index(res.distance) += 2*res.weight;
            }
        }
        return p_ax;
    };

    auto calc_wx = [&data_h, &data_x] (unsigned int imin, unsigned int imax) {
        container::Container1D<double> p_wx(constants::axes::d_axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) { // waters
            unsigned int j = 0;                      // exv
            for (; j+7 < data_x.get_size(); j+=8) {
                auto res = data_h[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3], data_x[j+4], data_x[j+5], data_x[j+6], data_x[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_wx.index(res.distance[k]) += res.weight[k];}
            }

            for (; j+3 < data_x.get_size(); j+=4) {
                auto res = data_h[i].evaluate(data_x[j], data_x[j+1], data_x[j+2], data_x[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_wx.index(res.distance[k]) += res.weight[k];}
            }

            for (; j < data_x.get_size(); ++j) {
                auto res = data_h[i].evaluate(data_x[j]);
                p_wx.index(res.distance) += res.weight;
            }
        }
        return p_wx;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    unsigned int job_size = settings::general::detail::job_size;
    BS::multi_future<container::Container1D<double>> xx;
    for (unsigned int i = 0; i < protein->atom_size(); i+=job_size) {
        xx.push_back(pool->submit(calc_xx, i, std::min(i+job_size, (unsigned int)protein->atom_size())));
    }
    BS::multi_future<container::Container1D<double>> ax;
    for (unsigned int i = 0; i < protein->water_size(); i+=job_size) {
        ax.push_back(pool->submit(calc_ax, i, std::min(i+job_size, (unsigned int)protein->water_size())));
    }
    BS::multi_future<container::Container1D<double>> wx;
    for (unsigned int i = 0; i < protein->water_size(); i+=job_size) {
        wx.push_back(pool->submit(calc_wx, i, std::min(i+job_size, (unsigned int)protein->water_size())));
    }

    //#################//
    // COLLECT RESULTS //
    //#################//
    auto p_xx_future = pool->submit(
        [&](){
            container::Container1D<double> p_xx(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : xx.get()) {
                std::transform(p_xx.begin(), p_xx.end(), tmp.begin(), p_xx.begin(), std::plus<double>());
            }
            return p_xx;
        }
    );

    auto p_ax_future = pool->submit(
        [&](){
            container::Container1D<double> p_ax(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hp.get()) {
                std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<double>());
            }
            return p_hp;
        }
    );

    auto p_wx_future = pool->submit(
        [&](){
            container::Container1D<double> p_hh(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hh.get()) {
                std::transform(p_hh.begin(), p_hh.end(), tmp.begin(), p_hh.begin(), std::plus<double>());
            }
            return p_hh;            
        }
    );

    pool.wait_for_tasks();
    auto p_pp = p_pp_future.get();
    auto p_hp = p_hp_future.get();
    auto p_hh = p_hh_future.get();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_pp[0] = std::accumulate(data_x_p.get_data_x().begin(), data_x_p.get_data_x().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );
    p_hh[0] = std::accumulate(data_x_h.get_data_x().begin(), data_x_h.get_data_x().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );

    // calculate p_tot    
    std::vector<double> p_tot(constants::axes::d_axis.bins, 0);
    for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot[i] = p_pp[i] + p_hh[i] + 2*p_hp[i];}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p_tot), Axis(0, max_bin*constants::axes::d_axis.width(), max_bin));
}