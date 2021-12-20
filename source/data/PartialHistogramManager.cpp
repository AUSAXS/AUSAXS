#include "data/Atom.h"
#include "data/Body.h"
#include "data/StateManager.h"
#include "ScatteringHistogram.h"
#include "data/PartialHistogramManager.h"

CompactCoordinates::CompactCoordinates(const Body& body) : size(body.protein_atoms.size()), data(4*size) {
    for (size_t i = 0; i < size; i++) {
        const Atom& a = body.protein_atoms[i]; 
        data[4*i] = a.coords.x;
        data[4*i+1] = a.coords.y;
        data[4*i+2] = a.coords.z;
        data[4*i+3] = a.effective_charge*a.occupancy;
    }
}

CompactCoordinates::CompactCoordinates(const vector<Hetatom>& atoms) : size(atoms.size()), data(4*size) {
    for (size_t i = 0; i < size; i++) {
        const Hetatom& a = atoms[i]; 
        data[4*i] = a.coords.x;
        data[4*i+1] = a.coords.y;
        data[4*i+2] = a.coords.z;
        data[4*i+3] = a.effective_charge*a.occupancy;
    }
}

PartialHistogramManager::PartialHistogramManager(vector<Body>& bodies, const vector<Hetatom>& hydration_atoms) 
    : size(bodies.size()), statemanager(size), coords_p(size), bodies(bodies), hydration_atoms(hydration_atoms), 
      partials_pp(size, vector<PartialHistogram>(size)), partials_hp(size) 
    {
        for (size_t i = 0; i < size; i++) {bodies[i].register_probe(statemanager.get_probe(i));}
    }

/**
 * @brief This initializes some necessary variables and precalculates the internal distances between atoms in each body.
 */
void PartialHistogramManager::initialize() {
    // generous sizes - 1000Å should be enough for just about any structure
    double width = setting::axes::scattering_intensity_plot_binned_width;
    vector<int> axes = {int(1000/width), 0, 1000}; 
    vector<double> p_base(axes[0], 0);

    for (size_t n = 0; n < size; n++) {
        // create more efficient access to the necessary variables
        CompactCoordinates current(bodies[n]);

        // calculate internal distances between atoms
        for (size_t i = 0; i < current.size; i++) {
            for (size_t j = i+1; j < current.size; j++) {
                float weight = current.data[4*i+3]*current.data[4*j+3];
                float dx = current.data[4*i] - current.data[4*j];
                float dy = current.data[4*i+1] - current.data[4*j+1];
                float dz = current.data[4*i+2] - current.data[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_base[dist/width] += 2*weight;
            }
        }

        // calculate self-correlation
        for (size_t i = 0; i < current.size; i++) {p_base[0] += current.data[4*i+3]*current.data[4*i+3];}

        // store the coordinates for later
        coords_p[n] = std::move(current);
    }
    master = MasterHistogram(p_base, axes);
}

ScatteringHistogram PartialHistogramManager::calculate_all() {
    Histogram total = calculate();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    vector<double> p_hh = partials_hh.p;
    
    vector<double> p_pp(total.axes[0]), p_hp(total.axes[0]);
    // iterate through all partial histograms in the upper triangle
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < i; j++) {
            PartialHistogram& current = partials_pp[i][j];

            // iterate through each entry in the partial histogram
            for (size_t k = 0; k < current.p.size(); k++) {
                p_pp[k] += current.p[k]; // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (size_t i = 0; i < size; i++) {
        PartialHistogram& current = partials_hp[i];

        // iterate through each entry in the partial histogram
        for (size_t k = 0; k < current.p.size(); k++) {
            p_hp[k] += current.p[k]; // add to p_pp
        }
    }

    return ScatteringHistogram(p_pp, p_hh, p_hp, total.p, total.axes);
}

Histogram PartialHistogramManager::calculate() {
    if (master.p.size() == 0) {initialize();} // check if this object has already been initialized

    // first we have to update the compact coordinate representations
    const vector<bool> modified_state = statemanager.get_modified_bodies();
    for (size_t i = 0; i < size; i++) {
        if (modified_state[i]) {
            coords_p[i] = CompactCoordinates(bodies[i]); // REMOVE - UNNECESSARY ON FIRST ITERATION
        }
    }

    // check if the hydration layer was modified
    if (statemanager.get_modified_hydration()) {
        coords_h = CompactCoordinates(hydration_atoms); // if so, first update the compact coordinate representation
        calc_hh(); // then update the partial histogram

        // iterate through the lower triangle
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < i; j++) {
                if (modified_state[i] || modified_state[j]) {
                    calc_pp(i, j);
                }
            }
            calc_hp(i); // we then update its partial histograms
        }
    }

    // if the hydration layer was not modified
    else {
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < i; j++) {
                if (modified_state[i] || modified_state[j]) { // if either of the two bodies were modified
                    calc_pp(i, j); // recalculate their partial histogram
                }
            }
            if (modified_state[i]) { // if a body was modified
                calc_hp(i); // update its partial histogram with the hydration layer
            }
        }
    }
    return Histogram(master.p, master.axes);
}

void PartialHistogramManager::calc_pp(const size_t& n, const size_t& m) {
    const double& width = setting::axes::scattering_intensity_plot_binned_width;
    const vector<int>& axes = master.axes; 

    CompactCoordinates& coords_n = coords_p[n];
    CompactCoordinates& coords_m = coords_p[m];
    vector<double> p_pp(axes[0], 0);
    for (size_t i = 0; i < coords_n.size; i++) {
        for (size_t j = 0; j < coords_m.size; j++) {
            float weight = coords_n.data[4*i+3]*coords_m.data[4*j+3];
            float dx = coords_n.data[4*i] - coords_m.data[4*j];
            float dy = coords_n.data[4*i+1] - coords_m.data[4*j+1];
            float dz = coords_n.data[4*i+2] - coords_m.data[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }
    master -= partials_pp[n][m];
    partials_pp[n][m].p = std::move(p_pp);
    master += partials_pp[n][m];
}


void PartialHistogramManager::calc_pp(const size_t& index) {
    const double& width = setting::axes::scattering_intensity_plot_binned_width;
    const vector<int>& axes = master.axes; 
    CompactCoordinates& coords_i = coords_p[index];

    // we do not want to calculate the self-correlation, so we have to skip entry 'index'
    for (size_t n = 0; n < index; n++) { // loop from (0, index]
        CompactCoordinates& coords_j = coords_p[n];
        vector<double> p_pp(axes[0], 0);
        for (size_t i = 0; i < coords_i.size; i++) {
            for (size_t j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[4*i+3]*coords_j.data[4*j+3];
                float dx = coords_i.data[4*i] - coords_j.data[4*j];
                float dy = coords_i.data[4*i+1] - coords_j.data[4*j+1];
                float dz = coords_i.data[4*i+2] - coords_j.data[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp[index][n];
        partials_pp[index][n].p = std::move(p_pp);
        master += partials_pp[index][n];
    }

    for (size_t n = index+1; n < size; n++) { // loop from (index, size]
        CompactCoordinates& coords_j = coords_p[n];
        vector<double> p_pp(axes[0], 0);
        for (size_t i = 0; i < coords_i.size; i++) {
            for (size_t j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[4*i+3]*coords_j.data[4*j+3];
                float dx = coords_i.data[4*i] - coords_j.data[4*j];
                float dy = coords_i.data[4*i+1] - coords_j.data[4*j+1];
                float dz = coords_i.data[4*i+2] - coords_j.data[4*j+2];
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp[index][n];
        partials_pp[index][n].p = std::move(p_pp);
        master += partials_pp[index][n];
    }
}

void PartialHistogramManager::calc_hp(const size_t& index) {
    const double& width = setting::axes::scattering_intensity_plot_binned_width;
    const vector<int>& axes = master.axes; 
    vector<double> p_hp(axes[0], 0);

    CompactCoordinates& coords = coords_p[index];
    for (size_t i = 0; i < coords.size; i++) {
        for (size_t j = 0; j < coords_h.size; j++) {
            float weight = coords.data[4*i+3]*coords_h.data[4*j+3];
            float dx = coords.data[4*i] - coords_h.data[4*j];
            float dy = coords.data[4*i+1] - coords_h.data[4*j+1];
            float dz = coords.data[4*i+2] - coords_h.data[4*j+2];
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
    const vector<int>& axes = master.axes; 
    vector<double> p_hh(axes[0], 0);

    // calculate internal distances for the hydration layer
    coords_h = CompactCoordinates(hydration_atoms);
    for (size_t i = 0; i < hydration_atoms.size(); i++) {
        for (size_t j = i+1; j < hydration_atoms.size(); j++) {
            float weight = coords_h.data[4*i+3]*coords_h.data[4*j+3];
            float dx = coords_h.data[4*i] - coords_h.data[4*j];
            float dy = coords_h.data[4*i+1] - coords_h.data[4*j+1];
            float dz = coords_h.data[4*i+2] - coords_h.data[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (size_t i = 0; i < hydration_atoms.size(); i++) {p_hh[0] += coords_h.data[4*i+3]*coords_h.data[4*i+3];}

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.p = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
}