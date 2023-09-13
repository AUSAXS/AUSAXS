#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <data/state/StateManager.h>
#include <hist/HistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <settings/HistogramSettings.h>

using namespace hist;

HistogramManager::HistogramManager(Protein* protein) : BodyTracker(protein), protein(protein) {}

HistogramManager::HistogramManager(const HistogramManager& hm) : BodyTracker(hm.protein), protein(hm.protein) {}

HistogramManager::~HistogramManager() = default;

std::unique_ptr<DistanceHistogram> HistogramManager::calculate() {return calculate_all();}

std::unique_ptr<CompositeDistanceHistogram> HistogramManager::calculate_all() {
    auto atoms = protein->get_atoms();
    auto waters = protein->get_waters();

    double width = settings::axes::distance_bin_width;
    Axis axes(0, settings::axes::max_distance, settings::axes::max_distance/width); 
    std::vector<double> p_pp(axes.bins, 0);
    std::vector<double> p_hh(axes.bins, 0);
    std::vector<double> p_hp(axes.bins, 0);
    std::vector<double> p_tot(axes.bins, 0);

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time
    std::vector<float> data_p(atoms.size()*4);
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        const Atom& a = atoms[i]; 
        data_p[4*i] = a.coords.x();
        data_p[4*i+1] = a.coords.y();
        data_p[4*i+2] = a.coords.z();
        data_p[4*i+3] = a.effective_charge*a.occupancy;
    }

    std::vector<float> data_h(waters.size()*4);
    for (unsigned int i = 0; i < waters.size(); ++i) {
        const Water& a = waters[i]; 
        data_h[4*i] = a.coords.x();
        data_h[4*i+1] = a.coords.y();
        data_h[4*i+2] = a.coords.z();
        data_h[4*i+3] = a.effective_charge*a.occupancy;
    }

    // calculate p-p distances
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = i+1; j < atoms.size(); ++j) {
            float weight = data_p[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_p[4*i] - data_p[4*j];
            float dy = data_p[4*i+1] - data_p[4*j+1];
            float dz = data_p[4*i+2] - data_p[4*j+2];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        p_pp[0] += data_p[4*i+3]*data_p[4*i+3];
    }

    for (unsigned int i = 0; i < waters.size(); ++i) {
        // calculate h-h distances
        for (unsigned int j = i+1; j < waters.size(); ++j) {
            float weight = data_h[4*i+3]*data_h[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_h[4*j];
            float dy = data_h[4*i+1] - data_h[4*j+1];
            float dz = data_h[4*i+2] - data_h[4*j+2];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }

        // calculate h-p distances
        for (unsigned int j = 0; j < atoms.size(); ++j) {
            float weight = data_h[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_p[4*j];
            float dy = data_h[4*i+1] - data_p[4*j+1];
            float dz = data_h[4*i+2] - data_p[4*j+2];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < waters.size(); ++i) {
        p_hh[0] += data_h[4*i+3]*data_h[4*i+3];
    }

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
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
    for (unsigned int i = 0; i < max_bin; ++i) {p_tot[i] = p_pp[i] + p_hh[i] + p_hp[i];}

    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hh), std::move(p_hp), std::move(p_tot), axes);
}