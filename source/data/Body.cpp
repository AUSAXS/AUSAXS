// includes
#include <vector>
#include <map>
#include <boost/format.hpp>
#include <utility>
#include <algorithm>

// my own includes
#include "data/Atom.h"
#include "hydrate/Grid.h"
#include "constants.h"
#include "data/Body.h"
#include "settings.h"
#include "math/Matrix.h"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::unique_ptr;

void Body::save(string path) {file->write(path);}

void Body::calc_histogram() {
    if (!updated_charge) {
        update_effective_charge(); // update the effective charge of all proteins. We have to do this since it affects the weights. 
    }

    // generous sizes - 1000Å should be enough for just about any structure
    double width = setting::axes::scattering_intensity_plot_binned_width;
    Axis axes = Axis(1000/width, 0, 1000); 
    vector<double> p_pp(axes.bins, 0);
    vector<double> p_hh(axes.bins, 0);
    vector<double> p_hp(axes.bins, 0);
    vector<double> p_tot(axes.bins, 0);

    // extremely wasteful to calculate this from scratch every time
    std::vector<float> data_p(protein_atoms.size()*4);
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        const Atom& a = protein_atoms[i]; 
        data_p[4*i] = a.coords.x;
        data_p[4*i+1] = a.coords.y;
        data_p[4*i+2] = a.coords.z;
        data_p[4*i+3] = a.effective_charge*a.occupancy;
    }

    std::vector<float> data_h(hydration_atoms.size()*4);
    for (size_t i = 0; i < hydration_atoms.size(); i++) {
        const Hetatom& a = hydration_atoms[i]; 
        data_h[4*i] = a.coords.x;
        data_h[4*i+1] = a.coords.y;
        data_h[4*i+2] = a.coords.z;
        data_h[4*i+3] = a.effective_charge*a.occupancy;
    }

    // calculate p-p distances
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        for (size_t j = i+1; j < protein_atoms.size(); j++) {
            float weight = data_p[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_p[4*i] - data_p[4*j];
            float dy = data_p[4*i+1] - data_p[4*j+1];
            float dz = data_p[4*i+2] - data_p[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (size_t i = 0; i < protein_atoms.size(); i++) {p_pp[0] += data_p[4*i+3]*data_p[4*i+3];}

    for (size_t i = 0; i < hydration_atoms.size(); i++) {
        // calculate h-h distances
        for (size_t j = i+1; j < hydration_atoms.size(); j++) {
            float weight = data_h[4*i+3]*data_h[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_h[4*j];
            float dy = data_h[4*i+1] - data_h[4*j+1];
            float dz = data_h[4*i+2] - data_h[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }

        // calculate h-p distances
        for (size_t j = 0; j < protein_atoms.size(); j++) {
            float weight = data_h[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_h[4*i] - data_p[4*j];
            float dy = data_h[4*i+1] - data_p[4*j+1];
            float dz = data_h[4*i+2] - data_p[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (size_t i = 0; i < hydration_atoms.size(); i++) {p_hh[0] += data_h[4*i+3]*data_h[4*i+3];}

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

    histogram = std::make_shared<ScatteringHistogram>(p_pp, p_hh, p_hp, p_tot, axes);
}

void Body::generate_new_hydration() {
    // delete the old hydration layer
    hydration_atoms = vector<Hetatom>();

    // move protein to center of mass
    center();

    // create the grid and hydrate it
    create_grid();
    hydration_atoms = grid->hydrate();
}

shared_ptr<Grid> Body::get_grid() {
    if (grid == nullptr) {create_grid();}
    return grid;
}

void Body::generate_volume_file(const string& path) {
    vector<vector<vector<char>>>& g = grid->grid;
    vector<Atom> filled;
    for (size_t i = 0; i < g.size(); i++) {
        for (size_t j = 0; j < g[0].size(); j++) {
            for (size_t k = 0; k < g[0][0].size(); k++) {
                if (g[i][j][k] != 0) {
                    Atom a(1, "CA", " ", "LEU", "A", 1, "", Vector3(i, j, k), 1, 0, "C", "");
                    filled.push_back(a);
                }
            }
        }
    }
    protein_atoms = filled;
    hydration_atoms = vector<Hetatom>();
    save(path);
    exit(0);
}

void Body::center() {
    if (!centered && setting::protein::center) {
        translate(-get_cm());
        centered = true;
    }
}

Vector3 Body::get_cm() const {
    Vector3 cm;
    double M = 0; // total mass
    auto weighted_sum = [&cm, &M] (auto& atoms) {
        for (auto const& a : atoms) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }
    };
    weighted_sum(protein_atoms);
    weighted_sum(hydration_atoms);
    return cm/M;
}

double Body::get_volume_acids() const {
    double v = 0;
    int cur_seq = 0; // sequence number of current acid
    for (auto const& a : protein_atoms) {
        int a_seq = a.resSeq; // sequence number of current atom
        if (cur_seq != a_seq) { // check if we are still dealing with the same acid
            cur_seq = a_seq; // if not, update our current sequence number
            v += constants::volume::get.at(a.resName); // and add its volume to the running total
        }
    }
    return v;
}

double Body::get_volume_grid() {
    if (grid == nullptr) {create_grid();}
    return grid->get_volume();
}

void Body::create_grid() {
    Axis3D axes(setting::grid::axes, setting::grid::width);
    grid = std::make_shared<Grid>(axes, setting::grid::width); 
    grid->add(protein_atoms);
    // grid->add(hydration_atoms);
}

shared_ptr<ScatteringHistogram> Body::get_histogram() {
    if (histogram == nullptr) {calc_histogram();}
    return histogram;
}

void Body::translate(const Vector3& v) {
    signal->state_change();

    auto move = [&v] (auto& atoms) {
        for (auto& a : atoms) {
            a.translate(v);
        }
    };
    move(protein_atoms);
    move(hydration_atoms);
}

void Body::rotate(const Matrix& R) {
    for (auto& atom : protein_atoms) {
        atom.coords.rotate(R);
    }

    for (auto& atom : hydration_atoms) {
        atom.coords.rotate(R);
    }
}

void Body::rotate(const double, const double, const double) {
    signal->state_change();
    cout << "Not implemented yet. " << endl;
    exit(1);
}

void Body::rotate(const Vector3& axis_arg, const double angle) {
    signal->state_change();

    // we use the Euler-Rodrigues formulation
    Vector3 axis = axis_arg.normalize_copy();
    double a = cos(angle/2);
    double b = sin(angle/2);
    double c = b;
    double d = b;
    b *= axis.x;
    c *= axis.y;
    d *= axis.z;

    double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

    Matrix R{{aa+bb-cc-dd, 2*(bc-ad),   2*(bd+ac)}, 
             {2*(bc+ad),   aa+cc-bb-dd, 2*(cd-ab)},
             {2*(bd-ac),   2*(cd+ab),   aa+dd-bb-cc}};

    rotate(R);
}

void Body::update_effective_charge(const double charge) {
    signal->state_change();
    std::for_each(protein_atoms.begin(), protein_atoms.end(), [&charge] (Atom& a) {a.add_effective_charge(charge);});
    updated_charge = true;
}

void Body::update_effective_charge() {
    signal->state_change();

    double displaced_vol = get_volume_grid();
    double displaced_charge = constants::charge::density::water*displaced_vol;
    // cout << "Displaced volume: " << displaced_vol << ", displaced charge: " << displaced_charge << endl;

    double charge_per_atom = -displaced_charge/protein_atoms.size();
    cout << "Added " << charge_per_atom << " additional charge to each protein atom (N: " << protein_atoms.size() << ")." << endl;

    std::for_each(protein_atoms.begin(), protein_atoms.end(), [&charge_per_atom] (Atom& a) {a.add_effective_charge(charge_per_atom);});
    updated_charge = true;
}

double Body::get_mass() const {
    double M = 0;
    std::for_each(protein_atoms.begin(), protein_atoms.end(), [&M] (const Atom& a) {M += a.get_mass();});
    std::for_each(hydration_atoms.begin(), hydration_atoms.end(), [&M] (const Hetatom& a) {M += a.get_mass();});
    return M;
}

Body& Body::operator=(const Body& rhs) {
    file = rhs.file;
    protein_atoms = file->protein_atoms;
    hydration_atoms = file->hydration_atoms;
    grid = rhs.grid;
    histogram = rhs.histogram;
    uid = rhs.uid;
    return *this;
}

bool Body::operator==(const Body& rhs) const {
    return uid == rhs.uid;
}