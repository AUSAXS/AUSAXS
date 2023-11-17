#pragma once

#include <data/Molecule.h>
#include <utility/Concepts.h>

/**
 * @brief Debug molecule that allows scaling the volume.
 */
class DebugMolecule : public data::Molecule {
    public:
        using Molecule::Molecule;
        double get_volume_grid() const override {return volume_scaling*data::Molecule::get_volume_grid();}
        void set_volume_scaling(double scaling) {volume_scaling = scaling;}

    private:
        double volume_scaling = 1;
};

template<container_type T1, container_type T2>
bool compare_hist(T1 p1, T2 p2) {
    int pmin = std::min<int>(p1.size(), p2.size());
    for (int i = 0; i < pmin; i++) {
        if (!utility::approx(p1[i], p2[i], 1e-6, 0)) {
            std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
            return false;
        }
    }

    return true;
}

template<typename T>
void set_unity_charge(T& protein) {
    // set the weights to 1 so we can analytically determine the result
    // waters
    for (auto& atom : protein.get_waters()) {
        atom.set_effective_charge(1);
    }
    // atoms
    for (auto& body : protein.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.set_effective_charge(1);
        }
    }
}

// calculation: 8 points
//          1 line  of length 0
//          3 lines of length 2
//          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
//          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
//
// calculation: 1 center point
//          1 line  of length 0
//          16 lines of length sqrt(3) = 1.73 (counting both directions)
//
// sum:
//          9 line  of length 0
//          16 lines of length sqrt(3)
//          24 lines of length 2
//          24 lines of length sqrt(8)
//          8 lines of length sqrt(12)
auto width = constants::axes::d_axis.width();
std::vector<double> d = {
    0, 
    constants::axes::d_vals[std::round(std::sqrt(3)/width)], 
    constants::axes::d_vals[std::round(2./width)], 
    constants::axes::d_vals[std::round(std::sqrt(8)/width)], 
    constants::axes::d_vals[std::round(std::sqrt(12)/width)]
};

std::vector<double> d_exact = {
    0, 
    std::sqrt(3), 
    2, 
    std::sqrt(8), 
    std::sqrt(12)
};