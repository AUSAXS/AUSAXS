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
            double sum1 = 0, sum2 = 0;
            for (int offset = -1; offset <= 1; ++offset) {
                int idx = i + offset;
                if (0 <= idx && idx < pmin) {
                    sum1 += p1[idx];
                    sum2 += p2[idx];    
                }
            }
            if (!utility::approx(sum1, sum2, abs, rel)) {
                std::cout << "Failed on index " << i << ". Window values: " << std::endl;
                for (int offset = -1; offset <= 1; ++offset) {
                    int idx = i + offset;
                    if (0 <= idx && idx < pmin) {
                        std::cout << "\t" << p1[idx] << ", " << p2[idx] << std::endl;
                    }
                }
                return false;
            }
        }
    }
    return true;
}

template<typename T>
void set_unity_charge(T& protein) {
    // set the weights to 1 so we can analytically determine the result
    for (auto& body : protein.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.weight() = 1;
        }
        auto w = body.get_waters();
        if (!w.has_value()) {continue;}
        for (auto& water : w->get()) {
            water.weight() = 1;
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

#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <catch2/catch_test_macros.hpp>

// This function uses template magic to invoke a given function template for all histogram manager variants, including weighted/unweighted and variable/fixed bin width versions.
// This is PART ONE of a two-part function, where this part handles the histogram managers that only vary in one template parameter.
template<typename F, typename ...Args>
void invoke_for_all_histogram_manager_variants_one(F&& f, Args&&... args) {
    SECTION("HistogramManagerMTFFGrid") {f.template operator()<hist::HistogramManagerMTFFGrid>((args)...);}
    SECTION("HistogramManagerMTFFGridSurface") {f.template operator()<hist::HistogramManagerMTFFGridSurface>((args)...);}
    SECTION("HistogramManagerMTFFGridScalableExv") {f.template operator()<hist::HistogramManagerMTFFGridScalableExv>((args)...);}
}

// This function uses template magic to invoke a given function template for all histogram manager variants, including weighted/unweighted and variable/fixed bin width versions.
// This is PART TWO of a two-part function, where this part handles the histogram managers that vary in two template parameters.
template<typename F, typename ...Args>
void invoke_for_all_histogram_manager_variants_two(F&& f, Args&&... args) {
    SECTION("HistogramManager") {f.template operator()<hist::HistogramManager>((args)...);}
    SECTION("HistogramManagerMT") {f.template operator()<hist::HistogramManagerMT>((args)...);}
    SECTION("HistogramManagerMTFFAvg") {f.template operator()<hist::HistogramManagerMTFFAvg>((args)...);}
    SECTION("HistogramManagerMTFFExplicit") {f.template operator()<hist::HistogramManagerMTFFExplicit>((args)...);}
    SECTION("SymmetryManagerMT") {f.template operator()<hist::SymmetryManagerMT>((args)...);}
    SECTION("PartialHistogramManager") {f.template operator()<hist::PartialHistogramManager>((args)...);}
    SECTION("PartialHistogramManagerMT") {f.template operator()<hist::PartialHistogramManagerMT>((args)...);}
    SECTION("PartialSymmetryManagerMT") {f.template operator()<hist::PartialSymmetryManagerMT>((args)...);}
}

template<typename F1, typename F2, typename ...Args>
void invoke_for_all_histogram_manager_variants(F1&& f1, F2&& f2, Args&&... args) {
    invoke_for_all_histogram_manager_variants_one(std::forward<F1>(f1), (args)...);
    invoke_for_all_histogram_manager_variants_two(std::forward<F2>(f2), (args)...);
}

template<typename F2, typename ...Args>
void invoke_for_all_nongrid_histogram_manager_variants(F2&& f2, Args&&... args) {
    invoke_for_all_histogram_manager_variants_two(std::forward<F2>(f2), (args)...);
}

template<typename F1, typename ...Args>
void invoke_for_all_grid_histogram_manager_variants(F1&& f1, Args&&... args) {
    invoke_for_all_histogram_manager_variants_one(std::forward<F1>(f1), (args)...);
}