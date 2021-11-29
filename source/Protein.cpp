// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
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
    if (path.find(".xml") != string::npos) { // .xml file
        print_err("Error in Protein::Protein: .xml input files are not supported.");
    } else if (path.find(".pdb") != string::npos) { // .pdb file
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

    // calculate the internal distances for the protein atoms
    int n_pp = 0; // index counter
    int n_hh = 0;
    int n_hp = 0;
    vector<double> d_pp(pow(protein_atoms.size(), 2));
    vector<double> w_pp(pow(protein_atoms.size(), 2)); 
    vector<double> d_hh(pow(hydration_atoms.size(), 2)); 
    vector<double> w_hh(pow(hydration_atoms.size(), 2)); 
    vector<double> d_hp(hydration_atoms.size()*protein_atoms.size());
    vector<double> w_hp(hydration_atoms.size()*protein_atoms.size()); 

    // calculate p-p distances
    for (int i = 0; i < protein_atoms.size(); i++) {
        for (int j = 0; j < protein_atoms.size(); j++) {
            d_pp[n_pp] = protein_atoms[i]->distance(protein_atoms[j]);
            w_pp[n_pp] = protein_atoms[i]->get_effective_charge()*protein_atoms[j]->get_effective_charge()
                *protein_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_pp++;
        }
    }

    for (int i = 0; i < hydration_atoms.size(); i++) {
        // calculate h-h distances
        for (int j = 0; j < hydration_atoms.size(); j++) {
            d_hh[n_hh] = hydration_atoms[i]->distance(hydration_atoms[j]);
            w_hh[n_hh] = hydration_atoms[i]->get_effective_charge()*hydration_atoms[j]->get_effective_charge()
                *hydration_atoms[i]->get_occupancy()*hydration_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_hh++;
        }
        // calculate h-p distances
        for (int j = 0; j < protein_atoms.size(); j++) {
            d_hp[n_hp] = hydration_atoms[i]->distance(protein_atoms[j]);
            w_hp[n_hp] = hydration_atoms[i]->get_effective_charge()*protein_atoms[j]->get_effective_charge()
                *hydration_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_hp++;
        }
    }
    this->distances = std::make_shared<Distances>(this, d_pp, d_hh, d_hp, w_pp, w_hh, w_hp);
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
    for (int i = 0; i < g.size(); i++) {
        for (int j = 0; j < g[0].size(); j++) {
            for (int k = 0; k < g[0][0].size(); k++) {
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