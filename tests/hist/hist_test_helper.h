#pragma once

#include <constants/ConstantsAxes.h>
#include <utility/Utility.h>
#include <data/Molecule.h>
#include <utility/Concepts.h>

#include <iostream>
#include <cmath>

using namespace ausaxs;

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

/**
 * @brief Check that the two containers are exactly identical. 
 */
template<container_type T1, container_type T2>
bool compare_hist(T1 p1, T2 p2, double abs = 1e-6, double rel = 1e-3) {
    int pmin = std::min<int>(p1.size(), p2.size());
    for (int i = 0; i < pmin; ++i) {
        if (!utility::approx(p1[i], p2[i], abs, rel)) {
            std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
            return false;
        }
    }
    return true;
}

/**
 * @brief Check that the two containers are approximately identical. 
 *        Variations across bin edges are allowed, meaning if a given bin is off by some amount, the following bin is checked for the difference. 
 */
template<container_type T1, container_type T2>
bool compare_hist_approx(T1 p1, T2 p2, double abs = 1e-6, double rel = 1e-3) {
    int pmin = std::min<int>(p1.size(), p2.size());
    for (int i = 0; i < pmin; ++i) {
        if (!utility::approx(p1[i], p2[i], abs, rel)) {
            if (i+1 < pmin) {
                auto diffi = std::abs(p1[i] - p2[i]);
                auto diffi1 = std::abs(p1[i+1] - p2[i+1]);
                if (!utility::approx(diffi, diffi1, abs, rel)) {
                    std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
                    std::cout << "Difference to next bin: " << diffi << ", " << diffi1 << std::endl;
                    return false;
                } else {
                    ++i;
                }
            } else {
                std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
                return false;
            }
        }
    }
    return true;
}

template<typename T>
void set_unity_charge(T& protein) {
    // set the weights to 1 so we can analytically determine the result
    // waters
    for (auto& atom : protein.get_waters()) {
        atom.weight() = 1;
    }
    // atoms
    for (auto& body : protein.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.weight() = 1;
        }
    }
}

struct SimpleCube {
    inline static std::vector<data::AtomFF> get_atoms() {
        return {
            data::AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), data::AtomFF({-1, 1, -1}, form_factor::form_factor_t::C),
            data::AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), data::AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C),
            data::AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), data::AtomFF({-1, 1,  1}, form_factor::form_factor_t::C),
            data::AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), data::AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)
        };
    }

    inline static std::vector<data::Water> get_waters() {
        return {
            data::Water({-1, -1, -1}), data::Water({-1, 1, -1}),
            data::Water({ 1, -1, -1}), data::Water({ 1, 1, -1}),
            data::Water({-1, -1,  1}), data::Water({-1, 1,  1}),
            data::Water({ 1, -1,  1}), data::Water({ 1, 1,  1})
        };
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
    inline static auto width = constants::axes::d_axis.width();
    inline static std::vector<double> d = {
        0, 
        constants::axes::d_vals[std::round(std::sqrt(3)/width)], 
        constants::axes::d_vals[std::round(2./width)], 
        constants::axes::d_vals[std::round(std::sqrt(8)/width)], 
        constants::axes::d_vals[std::round(std::sqrt(12)/width)]
    };

    inline static std::vector<double> d_exact = {
        0, 
        std::sqrt(3), 
        2, 
        std::sqrt(8), 
        std::sqrt(12)
    };

    inline static auto check_default = [] (const std::vector<double>& p) {
        if (p.back() < 2) {
            std::cout << "Failed on size: expected last index larger than 2Å, got: " << p.back() << std::endl;
            return false;
        }
        for (unsigned int i = 0; i < p.size(); ++i) {
            if (p[i] != constants::axes::d_vals[i]) {
                std::cout << "Failed on index " << i << ": expected: " << constants::axes::d_vals[i] << ", got: " << p[i] << std::endl;
                return false;
            }
        }
        return true;
    };

    inline static auto check_exact = [] (const std::vector<double>& p) {
        for (auto e : d_exact) {
            if (1e-6 < std::abs(p[std::round(e/constants::axes::d_axis.width())]-e)) {
                std::cout << "Failed on index " << std::round(e/constants::axes::d_axis.width()) << ": expected: " << e << ", got: " << p[std::round(e/constants::axes::d_axis.width())] << std::endl;
                return false;
            }
        }
        return true;
    };
};