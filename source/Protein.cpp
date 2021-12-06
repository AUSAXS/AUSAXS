// includes
#include <vector>
#include <map>
#include <boost/format.hpp>
#include <utility>
#include <algorithm>

// my own includes
#include "data/Atom.h"
#include "hydrate/Grid.h"
#include "data/PDB_file.cpp"
#include "data/constants.h"
#include "Protein.h"
#include "settings.h"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::unique_ptr;
using namespace ROOT;

Protein::Protein(string path) {
    // determine which kind of input file we're looking at
    if (path.find(".xml") != string::npos || path.find(".XML") != string::npos) { // .xml file
        print_err("Error in Protein::Protein: .xml input files are not supported.");
    } else if (path.find(".pdb") != string::npos || path.find(".PDB") != string::npos) { // .pdb file
        file = std::make_shared<PDB_file>(path);
    } else { // anything else - we cannot handle this
        print_err((format("Error in Protein::Protein: Invalid file extension of input file %1%.") % path).str());
        exit(1);
    }
    
    std::tie(protein_atoms, hydration_atoms) = file->get_atoms();
}

void Protein::save(string path) const {
    file->update(protein_atoms, hydration_atoms); // update the File backing this Protein with our new atoms
    file->write(path); // write to disk
}

void Protein::calc_distances() {
    update_effective_charge(); // update the effective charge of all proteins. We have to do this since it affects the weights. 

    // generous sizes - 1000Å should be enough for just about any structure
    vector<int> axes = {(int) (1000/setting::axes::scattering_intensity_plot_binned_width), 0, 1000}; 
    int max_bin = 0; // keep track of the highest distance stored, so everything else can be removed. 

    vector<double> p_pp(axes[0], 0);
    vector<double> p_hh(axes[0], 0);
    vector<double> p_hp(axes[0], 0);
    vector<double> p_tot(axes[0], 0);
    double width = (double) (axes[2]-axes[1])/axes[0]; // very important to cast this operation to a double - divison by two ints

    // calculate p-p distances
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        for (size_t j = 0; j < protein_atoms.size(); j++) {
            double dist = protein_atoms[i]->distance(protein_atoms[j]);
            double weight = protein_atoms[i]->get_effective_charge()*protein_atoms[j]->get_effective_charge()
                *protein_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2

            int index = std::round(dist/width);
            max_bin = std::max(index, max_bin);
            p_pp[index] += weight;
            p_tot[index] += weight;
        }
    }

    for (size_t i = 0; i < hydration_atoms.size(); i++) {
        // calculate h-h distances
        for (size_t j = 0; j < hydration_atoms.size(); j++) {
            double dist = hydration_atoms[i]->distance(hydration_atoms[j]);
            double weight = hydration_atoms[i]->get_effective_charge()*hydration_atoms[j]->get_effective_charge()
                *hydration_atoms[i]->get_occupancy()*hydration_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            int index = std::round(dist/width);
            max_bin = std::max(index, max_bin);
            p_hh[index] += weight;
            p_tot[index] += weight;
        }
        // calculate h-p distances
        for (size_t j = 0; j < protein_atoms.size(); j++) {
            double dist = hydration_atoms[i]->distance(protein_atoms[j]);
            double weight = hydration_atoms[i]->get_effective_charge()*protein_atoms[j]->get_effective_charge()
                *hydration_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            int index = std::round(dist/width);
            max_bin = std::max(index, max_bin);
            p_hp[index] += weight;
            p_tot[index] += weight;
        }
    }

    axes = {(int) (max_bin/setting::axes::scattering_intensity_plot_binned_width), 0, max_bin}; 
    cout << axes[0] << ", " << axes[1] << ", " << axes[2] << endl;
    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    this->distances = std::make_shared<Distances>(p_pp, p_hh, p_hp, p_tot, axes);
}

void Protein::generate_new_hydration() {
    // delete the old hydration layer
    hydration_atoms = vector<shared_ptr<Hetatom>>();

    // move protein to center of mass
    Vector3 cm = get_cm();
    translate(-cm);

    grid = std::make_shared<Grid>(setting::grid::base_point, setting::grid::width, setting::grid::bins/setting::grid::width); 
    grid->add(protein_atoms);
    hydration_atoms = grid->hydrate();
}

void Protein::generate_volume_file(string path) {
    vector<vector<vector<char>>>& g = grid->grid;
    vector<shared_ptr<Atom>> filled;
    for (size_t i = 0; i < g.size(); i++) {
        for (size_t j = 0; j < g[0].size(); j++) {
            for (size_t k = 0; k < g[0][0].size(); k++) {
                if (g[i][j][k] != 0) {
                    shared_ptr<Atom> a = std::make_shared<Atom>(0, "C", "", "C", "", 1, "", Vector3(i, j, k), 1, 0, "C", "");
                    filled.push_back(a);
                }
            }
        }
    }
    protein_atoms = filled;
    hydration_atoms = vector<shared_ptr<Hetatom>>();
    save(path);
    exit(0);
}

Vector3 Protein::get_cm() const {
    Vector3 cm;
    double M = 0; // total mass
    auto weighted_sum = [&cm, &M] (auto atoms) {
        for (auto const& a : *atoms) {
            double m = a->get_mass();
            M += m;
            double x = a->get_x()*m;
            double y = a->get_y()*m;
            double z = a->get_z()*m;
            cm += Vector3(x, y, z);
        }
        cm[0] = cm[0]/M;
        cm[1] = cm[1]/M;
        cm[2] = cm[2]/M;
    };
    weighted_sum(&protein_atoms);
    weighted_sum(&hydration_atoms);
    return cm;
}

double Protein::get_volume_acids() const {
    double v = 0;
    int cur_seq = 0; // sequence number of current acid
    for (auto const& a : protein_atoms) {
        int a_seq = a->get_resSeq(); // sequence number of current atom
        if (cur_seq != a_seq) { // check if we are still dealing with the same acid
            cur_seq = a_seq; // if not, update our current sequence number
            v += constants::volume::get.at(a->get_resName()); // and add its volume to the running total
        }
    }
    return v;
}

double Protein::get_volume_grid() {
    if (grid == nullptr) {create_grid();}
    return grid->get_volume();
}

void Protein::create_grid() {
    grid = std::make_shared<Grid>(setting::grid::base_point, setting::grid::width, setting::grid::bins/setting::grid::width); 
    grid->add(protein_atoms);
    grid->add(hydration_atoms);
}

shared_ptr<Distances> Protein::get_distances() {
    if (distances == nullptr) {calc_distances();}
    return distances;
}

void Protein::translate(const Vector3& v) {
    auto move = [&v] (auto atoms) {
        for (auto const& a : *atoms) {
            a->translate(v);
        }
    };
    move(&protein_atoms);
    move(&hydration_atoms);
}

void Protein::update_effective_charge() {
    if (grid == nullptr) {create_grid();}
    double displaced_vol = grid->get_volume();
    double displaced_charge = constants::charge::density::water*displaced_vol;
    double charge_per_atom = -displaced_charge/protein_atoms.size();
    cout << "Added " << charge_per_atom << " additional charge to each protein atom." << endl;
    std::for_each(protein_atoms.begin(), protein_atoms.end(), [&charge_per_atom] (const shared_ptr<Atom>& a) {a->add_effective_charge(charge_per_atom);});
}

double Protein::get_mass() const {
    double M = 0;
    std::for_each(protein_atoms.begin(), protein_atoms.end(), [&M] (const shared_ptr<Atom>& a) {M += a->get_mass();});
    std::for_each(hydration_atoms.begin(), hydration_atoms.end(), [&M] (const shared_ptr<Hetatom>& a) {M += a->get_mass();});
    cout << "Protein mass is " << M*constants::unit::gm << endl;
    return M*constants::unit::gm;
}