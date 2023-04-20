#include <data/Atom.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <data/StateManager.h>
#include <hist/ScatteringHistogram.h>
#include <hist/HistogramManager.h>
#include <hist/HistogramSettings.h>

using namespace hist;

HistogramManager::HistogramManager(Protein* protein) : BodyTracker(protein), protein(protein) {}

HistogramManager::~HistogramManager() = default;

Histogram HistogramManager::calculate() {return calculate_all().p;}

ScatteringHistogram HistogramManager::calculate_all() {
    auto atoms = protein->atoms();
    auto waters = protein->waters();

    double width = settings::axes::distance_bin_width;
    Axis axes = Axis(settings::axes::max_distance/width, 0, settings::axes::max_distance); 
    std::vector<double> p_pp(axes.bins, 0);
    std::vector<double> p_hh(axes.bins, 0);
    std::vector<double> p_hp(axes.bins, 0);
    std::vector<double> p_tot(axes.bins, 0);

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time
    std::vector<float> data_p(atoms.size()*4);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        const Atom& a = atoms[i]; 
        data_p[4*i] = a.coords.x();
        data_p[4*i+1] = a.coords.y();
        data_p[4*i+2] = a.coords.z();
        data_p[4*i+3] = a.effective_charge*a.occupancy;
    }

    std::vector<float> data_h(waters.size()*4);
    for (unsigned int i = 0; i < waters.size(); i++) {
        const Water& a = waters[i]; 
        data_h[4*i] = a.coords.x();
        data_h[4*i+1] = a.coords.y();
        data_h[4*i+2] = a.coords.z();
        data_h[4*i+3] = a.effective_charge*a.occupancy;
    }

    // calculate p-p distances
    for (unsigned int i = 0; i < atoms.size(); i++) {
        for (unsigned int j = i+1; j < atoms.size(); j++) {
            float weight = data_p[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_p[4*i] - data_p[4*j];
            float dy = data_p[4*i+1] - data_p[4*j+1];
            float dz = data_p[4*i+2] - data_p[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < atoms.size(); i++) {
        p_pp[0] += data_p[4*i+3]*data_p[4*i+3];
    }

    for (unsigned int i = 0; i < waters.size(); i++) {
        // calculate h-h distances
        for (unsigned int j = i+1; j < waters.size(); j++) {
            float weight = data_h[4*i+3]*data_h[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_h[4*j];
            float dy = data_h[4*i+1] - data_h[4*j+1];
            float dz = data_h[4*i+2] - data_h[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }

        // calculate h-p distances
        for (unsigned int j = 0; j < atoms.size(); j++) {
            float weight = data_h[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_p[4*j];
            float dy = data_h[4*i+1] - data_p[4*j+1];
            float dz = data_h[4*i+2] - data_p[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < waters.size(); i++) {
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
    axes = Axis{max_bin, 0, max_bin*width}; 

    // calculate p_tot    
    for (int i = 0; i < max_bin; i++) {p_tot[i] = p_pp[i] + p_hh[i] + p_hp[i];}

    return ScatteringHistogram(p_pp, p_hh, p_hp, p_tot, axes);
}