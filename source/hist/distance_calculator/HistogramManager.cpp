#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <data/state/StateManager.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <settings/HistogramSettings.h>

using namespace hist;

HistogramManager::HistogramManager(Protein* protein) : BodyTracker(protein), protein(protein) {}

HistogramManager::HistogramManager(const HistogramManager& hm) : BodyTracker(hm.protein), protein(hm.protein) {}

HistogramManager::~HistogramManager() = default;

std::unique_ptr<DistanceHistogram> HistogramManager::calculate() {return calculate_all();}

std::unique_ptr<CompositeDistanceHistogram> HistogramManager::calculate_all() {
    double width = settings::axes::distance_bin_width;
    Axis axes(0, settings::axes::max_distance, settings::axes::max_distance/width); 
    std::vector<double> p_pp(axes.bins, 0);
    std::vector<double> p_hh(axes.bins, 0);
    std::vector<double> p_hp(axes.bins, 0);
    std::vector<double> p_tot(axes.bins, 0);

    hist::detail::CompactCoordinates data_p(protein->get_bodies());
    hist::detail::CompactCoordinates data_h = hist::detail::CompactCoordinates(protein->get_waters());

    // calculate p-p distances
    for (unsigned int i = 0; i < data_p.get_size(); ++i) {
        for (unsigned int j = i+1; j < data_p.get_size(); ++j) {
            float weight = data_p[i].w*data_p[j].w;
            float dx = data_p[i].x - data_p[j].x;
            float dy = data_p[i].y - data_p[j].y;
            float dz = data_p[i].z - data_p[j].z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < data_p.get_size(); ++i) {
        p_pp[0] += std::pow(data_p[i].w, 2);
    }

    for (unsigned int i = 0; i < data_h.get_size(); ++i) {
        // calculate h-h distances
        for (unsigned int j = i+1; j < data_h.get_size(); ++j) {
            float weight = data_h[i].w*data_h[j].w;
            float dx = data_h[i].x - data_h[j].x;
            float dy = data_h[i].y - data_h[j].y;
            float dz = data_h[i].z - data_h[j].z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }

        // calculate h-p distances
        for (unsigned int j = 0; j < data_p.get_size(); ++j) {
            float weight = data_h[i].w*data_p[j].w;
            float dx = data_h[i].x - data_p[j].x;
            float dy = data_h[i].y - data_p[j].y;
            float dz = data_h[i].z - data_p[j].z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < data_h.get_size(); ++i) {
        p_hh[0] += std::pow(data_h[i].w, 2);
    }

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = axes.bins-1; i >= 10; i--) {
        if (p_pp[i] != 0 || p_hh[i] != 0 || p_hp[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    axes = Axis(0, max_bin*width, max_bin); 

    // calculate p_tot
    for (unsigned int i = 0; i < max_bin; ++i) {p_tot[i] = p_pp[i] + p_hh[i] + 2*p_hp[i];}

    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p_tot), axes);
}