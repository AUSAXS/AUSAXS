#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/MoleculeSettings.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

auto exact = [] (const data::Molecule& molecule, const std::vector<double>& qvals) {
    container::Container2D<double> distances(molecule.get_atoms().size(), molecule.get_atoms().size());
    auto atoms = molecule.get_atoms();
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = 0; j < atoms.size(); ++j) {
            distances(i, j) = atoms[i].coordinates().distance(atoms[j].coordinates());
        }
    }

    std::vector<double> I(qvals);
    for (unsigned int q = 0; q < qvals.size(); ++q) {
        double sum = 0;
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                double qd = qvals[q]*distances(i, j);
                if (qd < 1e-9) {
                    sum += atoms[i].weight()*atoms[j].weight();
                } else {
                    sum += std::sin(qd)/(qd)*atoms[i].weight()*atoms[j].weight();
                }
            }
        }
        I[q] = sum*std::exp(-qvals[q]*qvals[q]);
    }
    return I;
};

// Even the ExactDebyeCalculator uses the custom vectorized code for performance, so we should test that it agrees with the naïve implementation.
TEST_CASE("ExactDebyeCalculator: works") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    std::string file = GENERATE("2epe", "c60", "diamond");
    SECTION(file) {
        data::Molecule protein("tests/files/" + file + ".pdb");
        protein.clear_hydration();

        auto qaxis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
        auto I_exact = exact(protein, qaxis);
        auto Iq = hist::exact_debye_transform(protein, qaxis);

        REQUIRE(compare_hist(I_exact, Iq, 1e-6, 1e-6));
    }
}

// Test that the ExactDebyeCalculator agrees exactly with the analytical result for a simple system
TEST_CASE("ExactDebyeCalculator: agrees with analytical result") {
    auto d = SimpleCube::d_exact;
    static auto test_func = [&] (const auto& q_axis) {
        std::vector<double> Iq_exp;
        Iq_exp.resize(q_axis.size(), 0);
        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double dsum =
                8 +
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) +
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) +
                8 *std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            Iq_exp[q] += dsum*std::exp(-q_axis[q]*q_axis[q])*std::pow(form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C).I0(), 2);
        }
        return Iq_exp;
    };

    settings::molecule::implicit_hydrogens = false;
    auto protein = Molecule({Body{SimpleCube::get_atoms()}});

    SECTION("default q-axis") {
        auto Iq_exp = test_func(constants::axes::q_vals);
        auto Iq = hist::exact_debye_transform(protein, constants::axes::q_axis.as_vector());
        REQUIRE(compare_hist(Iq_exp, Iq, 1e-6, 1e-6));
    }

    SECTION("custom q-axis") {
        std::vector<double> q_axis(100);
        for (unsigned int i = 0; i < q_axis.size(); ++i) {
            q_axis[i] = (i+1)*0.1;
        }
        auto Iq_exp = test_func(q_axis);
        auto Iq = hist::exact_debye_transform(protein, q_axis);
        REQUIRE(compare_hist(Iq_exp, Iq, 1e-6, 1e-6));
    }
}