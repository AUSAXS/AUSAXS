#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "data/Atom.h"
#include "hydrate/Grid.h"
#include "data/PDB_file.cpp"
#include "data/properties.h"
#include "Protein.h"

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

void Protein::save(string path) {
    file->update(protein_atoms, hydration_atoms); // update the File backing this Protein with our new atoms
    file->write(path); // write to disk
}

Distances Protein::calc_distances() {
    // calculate the internal distances for the protein atoms
    int n_pp = 0; // index counter
    int n_hh = 0;
    int n_hp = 0;
    vector<double> d_pp(protein_atoms.size()*(protein_atoms.size() - 1)/2); // n(n-1)/2 entries
    vector<double> w_pp(protein_atoms.size()*(protein_atoms.size() - 1)/2); // corresponding weights
    vector<double> d_hh(hydration_atoms.size()*(hydration_atoms.size() - 1)/2); // m(m-1)/2 total entries
    vector<double> w_hh(hydration_atoms.size()*(hydration_atoms.size() - 1)/2); 
    vector<double> d_hp(hydration_atoms.size()*protein_atoms.size()); // n*m entries
    vector<double> w_hp(hydration_atoms.size()*protein_atoms.size()); 

    // calculate p-p distances
    for (int i = 0; i < protein_atoms.size(); i++) {
        for (int j = i+1; j < protein_atoms.size(); j++) {
            d_pp[n_pp] = protein_atoms[i]->distance(protein_atoms[j]);
            w_pp[n_pp] = property::charge::get.at(protein_atoms[i]->get_element())*property::charge::get.at(protein_atoms[j]->get_element())
                *protein_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_pp++;
        }
    }

    for (int i = 0; i < hydration_atoms.size(); i++) {
        // calculate h-h distances
        for (int j = i+1; j < hydration_atoms.size(); j++) {
            d_hh[n_hh] = hydration_atoms[i]->distance(hydration_atoms[j]);
            w_hh[n_hh] = property::charge::get.at(hydration_atoms[i]->get_element())*property::charge::get.at(hydration_atoms[j]->get_element())
                *hydration_atoms[i]->get_occupancy()*hydration_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_hh++;
        }
        // calculate h-p distances
        for (int j = 0; j < protein_atoms.size(); j++) {
            d_hp[n_hp] = hydration_atoms[i]->distance(protein_atoms[j]);
            w_hp[n_hp] = property::charge::get.at(hydration_atoms[i]->get_element())*property::charge::get.at(protein_atoms[j]->get_element())
                *hydration_atoms[i]->get_occupancy()*protein_atoms[j]->get_occupancy(); // Z1*Z2*w1*w2
            n_hp++;
        }
    }
    return Distances(d_pp, d_hh, d_hp, w_pp, w_hh, w_hp);
}

void Protein::generate_new_hydration(int reduce = 3, double width = 1) {
    // delete the old hydration layer
    hydration_atoms = vector<shared_ptr<Hetatom>>();

    // move protein to center of mass
    TVector3 cm = get_cm();
    translate(-cm);

    // generate the 3D grid
    Grid grid({-250, -250, -250}, width, 501/width); 
    grid.add(&protein_atoms);
    hydration_atoms = grid.hydrate(reduce);

    // double width = 10; // what width to use? 10 is too large, but with smaller values our grid becomes incredibly large
    // auto[corner, bins] = generate_grid(width); // corner is the lower corner of our grid, and bins the number of bins in each dimension
    // cout << format("bins: (%1%, %2%, %3%)") % bins[0] % bins[1] % bins[2] << endl;
}

std::pair<TVector3, vector<int>> Protein::generate_grid(const double width) {
    // determine the size of our grid
    TVector3 high = get_cm();
    TVector3 low = high;
    auto update = [&low, &high] (auto atoms) {
        for (auto const& a : *atoms) {
            // update minimum vector
            if (a->get_x() < low.X()) low.SetX(a->get_x());
            if (a->get_y() < low.Y()) low.SetY(a->get_y());
            if (a->get_z() < low.Z()) low.SetZ(a->get_z());

            // update maximum vector
            if (a->get_x() > high.X()) high.SetX(a->get_x());
            if (a->get_y() > high.Y()) high.SetY(a->get_y());
            if (a->get_z() > high.Z()) high.SetZ(a->get_z());
        }
    };
    update(&protein_atoms);
    update(&hydration_atoms);

    // calculate the number of bins in each dimension and initialize the grid and occupancy vectors
    vector<int> bins = {int((high.X() - low.X())/width), int((high.Y() - low.Y())/width), int((high.Z() - low.Z())/width)};
    return std::make_pair(low, bins);
}

TVector3 Protein::get_cm() {
    TVector3 cm;
    double M = 0; // total mass
    auto weighted_sum = [&cm, &M] (auto atoms) {
        for (auto const& a : *atoms) {
            double m = a->get_atomic_weight();
            M += m;
            double x = a->get_x()*m;
            double y = a->get_y()*m;
            double z = a->get_z()*m;
            cm += TVector3(x, y, z);
        }
        cm[0] = cm[0]/M;
        cm[1] = cm[1]/M;
        cm[2] = cm[2]/M;
    };
    weighted_sum(&protein_atoms);
    weighted_sum(&hydration_atoms);
    return cm;
}

double Protein::get_volume() {
    double v = 0;
    int cur_seq = 0; // sequence number of current acid
    for (auto const& a : protein_atoms) {
        int a_seq = a->get_resSeq(); // sequence number of current atom
        if (cur_seq != a_seq) { // check if we are still dealing with the same acid
            cur_seq = a_seq; // if not, update our current sequence number
            v += property::volume::get.at(a->get_resName()); // and add its volume to the running total
        }
    }
    return v;
}

void Protein::translate(const TVector3 v) {
    auto move = [&v] (auto atoms) {
        for (auto const& a : *atoms) {
            a->translate(v);
        }
    };
    move(&protein_atoms);
    move(&hydration_atoms);
}