#include <data/Atom.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <data/StateManager.h>
#include <histogram/ScatteringHistogram.h>
#include <histogram/PartialHistogramManager.h>

using namespace hist;

using std::vector;

CompactCoordinates::CompactCoordinates(const Body& body) : size(body.get_protein_atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; i++) {
        const Atom& a = body.protein_atom(i); 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

CompactCoordinates::CompactCoordinates(const vector<Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; i++) {
        const Water& a = atoms[i]; 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

MasterHistogram::MasterHistogram(const vector<double>& p_base, const Axis& axis) : Histogram(p_base, axis), base(std::move(p_base)) {}

MasterHistogram& MasterHistogram::operator+=(const PartialHistogram& rhs) {
    p += rhs.p;
    return *this;
}

MasterHistogram& MasterHistogram::operator-=(const PartialHistogram& rhs) {
    p -= rhs.p;
    return *this;
}

PartialHistogramManager::PartialHistogramManager(Protein* protein) 
    : size(protein->bodies.size()), statemanager(size), coords_p(size), protein(protein), partials_pp(size, vector<PartialHistogram>(size)), partials_hp(size) 
    {
        for (unsigned int i = 0; i < size; i++) {protein->bodies[i].register_probe(statemanager.get_probe(i));}
    }

std::shared_ptr<StateManager::BoundSignaller> PartialHistogramManager::get_probe(unsigned int i) {return statemanager.get_probe(i);}

void PartialHistogramManager::signal_modified_hydration_layer() {statemanager.modified_hydration_layer();}

const StateManager& PartialHistogramManager::get_state_manager() const {return statemanager;}

StateManager& PartialHistogramManager::get_state_manager() {return statemanager;}

void PartialHistogramManager::calc_self_correlation(unsigned int index) {
    double width = setting::axes::scattering_intensity_plot_binned_width;
    CompactCoordinates current(protein->bodies[index]);

    // calculate internal distances between atoms
    vector<double> p_pp(master.axis.bins, 0);
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
    // generous sizes - 1000Å should be enough for just about any structure
    double width = setting::axes::scattering_intensity_plot_binned_width;
    Axis axis = Axis(1000/width, 0, 1000); 
    vector<double> p_base(axis.bins, 0);
    master = MasterHistogram(p_base, axis);

    partials_hh = PartialHistogram(axis);
    for (unsigned int n = 0; n < size; n++) {
        partials_hp[n] = PartialHistogram(axis);
        partials_pp[n][n] = PartialHistogram(axis);
        calc_self_correlation(n);

        for (unsigned int k = 0; k < n; k++) {
            partials_pp[n][k] = PartialHistogram(axis);
        }
    }
}

ScatteringHistogram PartialHistogramManager::calculate_all() {
    Histogram total = calculate();
    total.shorten_axis();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    vector<double> p_hh = partials_hh.p;
    vector<double> p_pp = master.base.p;
    vector<double> p_hp(total.axis.bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < i; j++) {
            PartialHistogram& current = partials_pp[i][j];

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < total.axis.bins; k++) {
                p_pp[k] += current.p[k]; // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < size; i++) {
        PartialHistogram& current = partials_hp[i];

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
    const vector<bool> externally_modified = statemanager.get_externally_modified_bodies();
    const vector<bool> internally_modified = statemanager.get_internally_modified_bodies();

    // check if the object has already been initialized
    if (__builtin_expect(master.p.size() == 0, false)) {
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
                coords_p[i] = CompactCoordinates(protein->bodies[i]);
            }
        }
    }

    // check if the hydration layer was modified
    if (statemanager.get_modified_hydration()) {
        coords_h = CompactCoordinates(protein->hydration_atoms); // if so, first update the compact coordinate representation
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

    CompactCoordinates& coords_n = coords_p[n];
    CompactCoordinates& coords_m = coords_p[m];
    vector<double> p_pp(master.axis.bins, 0);
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
    CompactCoordinates& coords_i = coords_p[index];

    // we do not want to calculate the self-correlation, so we have to skip entry 'index'
    for (size_t n = 0; n < index; n++) { // loop from (0, index]
        CompactCoordinates& coords_j = coords_p[n];
        vector<double> p_pp(master.axis.bins, 0);
        for (size_t i = 0; i < coords_i.size; i++) {
            for (size_t j = 0; j < coords_j.size; j++) {
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

    for (size_t n = index+1; n < size; n++) { // loop from (index, size]
        CompactCoordinates& coords_j = coords_p[n];
        vector<double> p_pp(master.axis.bins, 0);
        for (size_t i = 0; i < coords_i.size; i++) {
            for (size_t j = 0; j < coords_j.size; j++) {
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
    vector<double> p_hp(master.axis.bins, 0);

    CompactCoordinates& coords = coords_p[index];
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
    vector<double> p_hh(master.axis.bins, 0);

    // calculate internal distances for the hydration layer
    coords_h = CompactCoordinates(protein->hydration_atoms);
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