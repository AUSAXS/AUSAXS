#include <data/Atom.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <data/StateManager.h>
#include <hist/ScatteringHistogram.h>
#include <hist/PartialHistogramManager.h>

using namespace hist;

PartialHistogramManager::PartialHistogramManager(Protein* protein) 
    : size(protein->bodies.size()), statemanager(size), coords_p(size), protein(protein), partials_pp(size, std::vector<detail::PartialHistogram>(size)), partials_hp(size) 
    {
        for (unsigned int i = 0; i < size; i++) {protein->bodies[i].register_probe(statemanager.get_probe(i));}
    }

std::shared_ptr<StateManager::BoundSignaller> PartialHistogramManager::get_probe(unsigned int i) {return statemanager.get_probe(i);}

void PartialHistogramManager::signal_modified_hydration_layer() {statemanager.modified_hydration_layer();}

const StateManager& PartialHistogramManager::get_state_manager() const {return statemanager;}

StateManager& PartialHistogramManager::get_state_manager() {return statemanager;}

void PartialHistogramManager::calc_self_correlation(unsigned int index) {
    double width = setting::axes::scattering_intensity_plot_binned_width;
    detail::CompactCoordinates current(protein->bodies[index]);

    // calculate internal distances between atoms
    std::vector<double> p_pp(master.axis.bins, 0);
    for (unsigned int i = 0; i < current.size; i++) {
        for (unsigned int j = i+1; j < current.size; j++) {
            float weight = current.data[i].w*current.data[j].w;
            float dx = current.data[i].x - current.data[j].x;
            float dy = current.data[i].y - current.data[j].y;
            float dz = current.data[i].z - current.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < current.size; i++) {p_pp[0] += current.data[i].w*current.data[i].w;}

    // store the coordinates for later
    coords_p[index] = std::move(current);

    master.base -= partials_pp[index][index];
    master -= partials_pp[index][index];
    partials_pp[index][index].p = std::move(p_pp);
    master += partials_pp[index][index];
    master.base += partials_pp[index][index];
}

/**
 * @brief This initializes some necessary variables and precalculates the internal distances between atoms in each body.
 */
void PartialHistogramManager::initialize() {
    double width = setting::axes::scattering_intensity_plot_binned_width;
    Axis axis = Axis(setting::axes::max_distance/width, 0, setting::axes::max_distance); 
    std::vector<double> p_base(axis.bins, 0);
    master = detail::MasterHistogram(p_base, axis);

    partials_hh = detail::PartialHistogram(axis);
    for (unsigned int n = 0; n < size; n++) {
        partials_hp[n] = detail::PartialHistogram(axis);
        partials_pp[n][n] = detail::PartialHistogram(axis);
        calc_self_correlation(n);

        for (unsigned int k = 0; k < n; k++) {
            partials_pp[n][k] = detail::PartialHistogram(axis);
        }
    }
}

ScatteringHistogram PartialHistogramManager::calculate_all() {
    Histogram total = calculate();
    total.shorten_axis();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    std::vector<double> p_hh = partials_hh.p;
    std::vector<double> p_pp = master.base.p;
    std::vector<double> p_hp(total.axis.bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < i; j++) {
            detail::PartialHistogram& current = partials_pp[i][j];

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < total.axis.bins; k++) {
                p_pp[k] += current.p[k]; // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < size; i++) {
        detail::PartialHistogram& current = partials_hp[i];

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < total.axis.bins; k++) {
            p_hp[k] += current.p[k]; // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(total.axis.bins);
    p_pp.resize(total.axis.bins);

    return ScatteringHistogram(p_pp, p_hh, p_hp, std::move(total.p), total.axis);
}

Histogram PartialHistogramManager::calculate() {
    const std::vector<bool> externally_modified = statemanager.get_externally_modified_bodies();
    const std::vector<bool> internally_modified = statemanager.get_internally_modified_bodies();

    // check if the object has already been initialized
    if (master.p.size() == 0) [[unlikely]] {
        initialize(); 
    } 
    
    // if not, we must first check if the coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < size; i++) {
            if (internally_modified[i]) {
                // if the internal state was modified, we have to recalculate the self-correlation
                calc_self_correlation(i);
            } else if (externally_modified[i]) {
                // if the external state was modified, we have to update the coordinate representations
                coords_p[i] = detail::CompactCoordinates(protein->bodies[i]);
            }
        }
    }

    // check if the hydration layer was modified
    if (statemanager.get_modified_hydration()) {
        coords_h = detail::CompactCoordinates(protein->hydration_atoms); // if so, first update the compact coordinate representation
        calc_hh(); // then update the partial histogram

        // iterate through the lower triangle
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) {
                    calc_pp(i, j);
                }
            }
            calc_hp(i); // we then update its partial histograms
        }
    }

    // if the hydration layer was not modified
    else {
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) { // if either of the two bodies were modified
                    calc_pp(i, j); // recalculate their partial histogram
                }
            }
            if (externally_modified[i]) { // if a body was modified
                calc_hp(i); // update its partial histogram with the hydration layer
            }
        }
    }
    statemanager.reset();
    return Histogram(master.p, master.axis);
}

void PartialHistogramManager::calc_pp(unsigned int n, unsigned int m) {
    double width = setting::axes::scattering_intensity_plot_binned_width;

    detail::CompactCoordinates& coords_n = coords_p[n];
    detail::CompactCoordinates& coords_m = coords_p[m];
    std::vector<double> p_pp(master.axis.bins, 0);
    for (unsigned int i = 0; i < coords_n.size; i++) {
        for (unsigned int j = 0; j < coords_m.size; j++) {
            float weight = coords_n.data[i].w*coords_m.data[j].w;
            float dx = coords_n.data[i].x - coords_m.data[j].x;
            float dy = coords_n.data[i].y - coords_m.data[j].y;
            float dz = coords_n.data[i].z - coords_m.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }
    master -= partials_pp[n][m];
    partials_pp[n][m].p = std::move(p_pp);
    master += partials_pp[n][m];
}

void PartialHistogramManager::calc_pp(unsigned int index) {
    double width = setting::axes::scattering_intensity_plot_binned_width;
    detail::CompactCoordinates& coords_i = coords_p[index];

    // we do not want to calculate the self-correlation, so we have to skip entry 'index'
    for (unsigned int n = 0; n < index; n++) { // loop from (0, index]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = 0; i < coords_i.size; i++) {
            for (unsigned int j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[i].w*coords_j.data[j].w;
                float dx = coords_i.data[i].x - coords_j.data[j].x;
                float dy = coords_i.data[i].y - coords_j.data[j].y;
                float dz = coords_i.data[i].z - coords_j.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp[index][n];
        partials_pp[index][n].p = std::move(p_pp);
        master += partials_pp[index][n];
    }

    for (unsigned int n = index+1; n < size; n++) { // loop from (index, size]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = 0; i < coords_i.size; i++) {
            for (unsigned int j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[i].w*coords_j.data[j].w;
                float dx = coords_i.data[i].x - coords_j.data[j].x;
                float dy = coords_i.data[i].y - coords_j.data[j].y;
                float dz = coords_i.data[i].z - coords_j.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp[index][n];
        partials_pp[index][n].p = std::move(p_pp);
        master += partials_pp[index][n];
    }
}

void PartialHistogramManager::calc_hp(unsigned int index) {
    double width = setting::axes::scattering_intensity_plot_binned_width;
    std::vector<double> p_hp(master.axis.bins, 0);

    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.size; i++) {
        for (unsigned int j = 0; j < coords_h.size; j++) {
            float weight = coords.data[i].w*coords_h.data[j].w;
            float dx = coords.data[i].x - coords_h.data[j].x;
            float dy = coords.data[i].y - coords_h.data[j].y;
            float dz = coords.data[i].z - coords_h.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += 2*weight;
        }
    }

    master -= partials_hp[index]; // subtract the previous hydration histogram
    partials_hp[index].p = std::move(p_hp);
    master += partials_hp[index]; // add the new hydration histogram
}

void PartialHistogramManager::calc_hh() {
    const double& width = setting::axes::scattering_intensity_plot_binned_width;
    std::vector<double> p_hh(master.axis.bins, 0);

    // calculate internal distances for the hydration layer
    coords_h = detail::CompactCoordinates(protein->hydration_atoms);
    for (unsigned int i = 0; i < 4*protein->hydration_atoms.size(); i++) {
        for (unsigned int j = i+1; j < protein->hydration_atoms.size(); j++) {
            float weight = coords_h.data[i].w*coords_h.data[j].w;
            float dx = coords_h.data[i].x - coords_h.data[j].x;
            float dy = coords_h.data[i].y - coords_h.data[j].y;
            float dz = coords_h.data[i].z - coords_h.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < protein->hydration_atoms.size(); i++) {p_hh[0] += coords_h.data[i].w*coords_h.data[i].w;}

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.p = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
}

ScatteringHistogram PartialHistogramManager::calculate_slow() const {
    auto atoms = protein->atoms();
    auto waters = protein->waters();

    double width = setting::axes::scattering_intensity_plot_binned_width;
    Axis axes = Axis(setting::axes::max_distance/width, 0, setting::axes::max_distance); 
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