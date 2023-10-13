#include <hist/distance_calculator/HistogramManagerMTFF.h>
#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/FormFactor.h>
#include <data/Protein.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <utility/Container3D.h>
#include <utility/Container2D.h>
#include <utility/Container1D.h>

#include <BS_thread_pool.hpp>

using namespace hist;

HistogramManagerMTFF::HistogramManagerMTFF(HistogramManager& hm) : HistogramManager(hm) {}

HistogramManagerMTFF::~HistogramManagerMTFF() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFF::calculate() {return calculate_all();}

std::unique_ptr<CompositeDistanceHistogram> HistogramManagerMTFF::calculate_all() {
    double width = settings::axes::distance_bin_width;
    double Z_exv_avg = protein->get_excluded_volume()*constants::charge::density::water/protein->atom_size();
    std::cout << "Z_exv_avg = " << Z_exv_avg << std::endl;
    double Z_exv_avg2 = Z_exv_avg*Z_exv_avg;
    unsigned int excluded_volume_bin = static_cast<unsigned int>(hist::detail::form_factor_t::EXCLUDED_VOLUME);
    Axis axes(0, settings::axes::max_distance, settings::axes::max_distance/width); 

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time (class is not meant for serial use anyway?)
    hist::detail::CompactCoordinatesFF data_p(protein->get_bodies());
    hist::detail::CompactCoordinatesFF data_h = hist::detail::CompactCoordinatesFF(protein->get_waters());

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    BS::thread_pool pool(settings::general::threads);
    auto calc_pp = [&data_p, &axes, &width, &Z_exv_avg, &Z_exv_avg2, &excluded_volume_bin] (unsigned int imin, unsigned int imax) {
        Container3D<double> p_pp(detail::FormFactor::get_count(), detail::FormFactor::get_count(), axes.bins, 0); // ff_type1, ff_type2, distance
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = 0; j < data_p.get_size(); ++j) {
                float dx = data_p[i].x - data_p[j].x;
                float dy = data_p[i].y - data_p[j].y;
                float dz = data_p[i].z - data_p[j].z;
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                p_pp.index(data_p[i].ff_type, data_p[j].ff_type, dist/width) += data_p[i].w*data_p[j].w;
                p_pp.index(data_p[i].ff_type, excluded_volume_bin, dist/width) += data_p[i].w*Z_exv_avg;
                p_pp.index(excluded_volume_bin, excluded_volume_bin, dist/width) += Z_exv_avg2;
            }
        }
        return p_pp;
    };

    auto calc_hp = [&data_h, &data_p, &axes, &width, &Z_exv_avg, &excluded_volume_bin] (unsigned int imin, unsigned int imax) {
        Container2D<double> p_hp(detail::FormFactor::get_count(), axes.bins, 0); // ff_type, distance
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = 0; j < data_p.get_size(); ++j) {
                float dx = data_h[i].x - data_p[j].x;
                float dy = data_h[i].y - data_p[j].y;
                float dz = data_h[i].z - data_p[j].z;
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                p_hp.index(data_p[j].ff_type, dist/width) += data_h[i].w*data_p[j].w;
                p_hp.index(excluded_volume_bin, dist/width) += data_h[i].w*Z_exv_avg;
            }
        }
        return p_hp;
    };

    auto calc_hh = [&data_h, &axes, &width] (unsigned int imin, unsigned int imax) {
        Container1D<double> p_hh(axes.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = 0; j < data_h.get_size(); ++j) {
                float weight = data_h[i].w*data_h[j].w;
                float dx = data_h[i].x - data_h[j].x;
                float dy = data_h[i].y - data_h[j].y;
                float dz = data_h[i].z - data_h[j].z;
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                p_hh.index(dist/width) += weight;
            }
        }
        return p_hh;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    unsigned int job_size = settings::general::detail::job_size;
    BS::multi_future<Container3D<double>> pp;
    for (unsigned int i = 0; i < data_p.get_size(); i+=job_size) {
        pp.push_back(pool.submit(calc_pp, i, std::min(i+job_size, data_p.get_size())));
    }
    BS::multi_future<Container2D<double>> hp;
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
            Container3D<double> p_pp(detail::FormFactor::get_count(), detail::FormFactor::get_count(), axes.bins, 0); // ff_type1, ff_type2, distance
            for (const auto& tmp : pp.get()) {
                std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<double>());
            }
            return p_pp;
        }
    );

    auto p_hp_future = pool.submit(
        [&]() {
            Container2D<double> p_hp(detail::FormFactor::get_count(), axes.bins, 0); // ff_type, distance
            for (const auto& tmp : hp.get()) {
                std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<double>());
            }
            return p_hp;
        }
    );

    auto p_hh_future = pool.submit(
        [&]() {
            Container1D<double> p_hh(axes.bins, 0);
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
    // for (unsigned int i = 0; i < data_p.get_size(); ++i) {
    //     p_pp.index(data_p[i].ff_type, data_p[i].ff_type, 0) -= std::pow(data_p[i].w, 2);
    //     p_pp.index(excluded_volume_bin, excluded_volume_bin, 0) -= std::pow(Z_exv_avg, 2);
    // }
    // for (unsigned int i = 0; i < data_h.get_size(); ++i) {
    //     p_hh.index(0) -= std::pow(data_h[i].w, 2);
    // }

    std::vector<double> p_tot(axes.bins, 0);
    for (unsigned int i = 0; i < axes.bins; ++i) {
        for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count_without_excluded_volume(); ++ff2) {
                p_tot[i] += p_pp.index(ff1, ff2, i);
            }
            p_tot[i] += 2*p_hp.index(ff1, i);
        }
        p_tot[i] += p_hh.index(i);
    }

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (unsigned int i = axes.bins-1; i >= 10; --i) {
        if (p_tot[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    Container3D<double> p_pp_short(detail::FormFactor::get_count(), detail::FormFactor::get_count(), max_bin);
    Container2D<double> p_hp_short(detail::FormFactor::get_count(), max_bin);
    Container1D<double> p_hh_short(max_bin);
    for (unsigned int i = 0; i < max_bin; ++i) {
        for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
                p_pp_short.index(ff1, ff2, i) = p_pp.index(ff1, ff2, i);
            }
            p_hp_short.index(ff1, i) = p_hp.index(ff1, i);
        }
        p_hh_short.index(i) = p_hh.index(i);
    }
    p_tot.resize(max_bin);

    axes = Axis(0, max_bin*width, max_bin); 
    return std::make_unique<CompositeDistanceHistogramFF>(std::move(p_pp_short), std::move(p_hp_short), std::move(p_hh_short), std::move(p_tot), axes);
}
