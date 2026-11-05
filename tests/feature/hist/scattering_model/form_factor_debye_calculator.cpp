#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/ConstantsAxes.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <form_factor/neutron/NeutronFormFactor.h>
#include <form_factor/xray/XrayFormFactor.h>
#include <hist/scattering_model/FormFactorDebyeCalculator.h>
#include <hist/scattering_model/IScatteringModelExv.h>
#include <math/ConstexprMath.h>
#include <settings/ExvSettings.h>
#include <settings/FormFactorSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/MoleculeSettings.h>
#include <settings/ScatteringSettings.h>

#include <utility/Exceptions.h>

#include <cmath>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using ff_t = form_factor::form_factor_t;
using Radiation = settings::scattering::Radiation;

namespace {
    struct Scatterer {
        Vector3<double> r;
        ff_t type;
    };

    template<typename F>
    double evaluate(ff_t type, double q, F&& eval) {
        if (settings::scattering::radiation == Radiation::Neutron) {return eval(form_factor::neutron::get(type), q);}
        return eval(form_factor::xray::raw::get(type), q);
    }

    // the direct sum over all pairs of scatterers, with the form factors taken straight from the tables
    std::vector<double> direct_sum(const std::vector<Scatterer>& s, const std::vector<double>& q) {
        auto f = [] (ff_t t, double q) {return evaluate(t, q, [] (const auto& ff, double q) {return ff.evaluate(q);});};
        auto f_self = [] (ff_t t, double q) {return evaluate(t, q, [] (const auto& ff, double q) {return ff.evaluate_self(q);});};
        std::vector<double> I(q.size(), 0);
        for (unsigned int k = 0; k < q.size(); ++k) {
            for (unsigned int i = 0; i < s.size(); ++i) {
                I[k] += f_self(s[i].type, q[k]);
                for (unsigned int j = 0; j < s.size(); ++j) {
                    if (i == j) {continue;}
                    I[k] += f(s[i].type, q[k])*f(s[j].type, q[k])*math::sinc(q[k]*s[i].r.distance(s[j].r));
                }
            }
        }
        return I;
    }

    // the scale of the intensity at q, (sum_i |f_i|)^2, against which the deviations are measured since the intensity crosses zero
    std::vector<double> scale(const std::vector<Scatterer>& s, const std::vector<double>& q) {
        std::vector<double> scale(q.size(), 0);
        for (unsigned int k = 0; k < q.size(); ++k) {
            for (const auto& si : s) {scale[k] += std::abs(evaluate(si.type, q[k], [] (const auto& ff, double q) {return ff.evaluate(q);}));}
            scale[k] *= scale[k];
        }
        return scale;
    }

    // scatterers on a line at the marks of a Golomb ruler, so every distance is distinct and has a bin of its own. Weighted bins then
    // hold the exact distances, and the binned sum must match the direct one. The last two are waters.
    std::vector<Scatterer> golomb_scatterers() {
        constexpr std::array marks = {0, 1, 5, 12, 25, 27, 35, 41, 44};
        constexpr std::array types = {ff_t::C, ff_t::N, ff_t::O, ff_t::S, ff_t::Fe, ff_t::Zn, ff_t::I, ff_t::WATER, ff_t::WATER};
        std::vector<Scatterer> s;
        for (unsigned int i = 0; i < marks.size(); ++i) {s.push_back({{1.13*marks[i], 0.2*i, 0}, types[i]});}
        return s;
    }

    Molecule as_molecule(const std::vector<Scatterer>& s) {
        std::vector<AtomFF> atoms;
        std::vector<Water> waters;
        for (const auto& si : s) {
            if (si.type == ff_t::WATER) {waters.emplace_back(si.r);}
            else {atoms.emplace_back(si.r, si.type);}
        }
        return Molecule({Body{atoms, waters}});
    }

    // restores the global settings these tests change
    struct SettingsGuard {
        Radiation radiation = settings::scattering::radiation;
        settings::exv::ExvMethod exv = settings::exv::exv_method;
        int max_types = settings::form_factor::max_types;
        double min_fraction = settings::form_factor::min_fraction;
        settings::hist::WeightedBins weighted_bins = settings::hist::weighted_bins;
        ~SettingsGuard() {
            settings::scattering::radiation = radiation;
            settings::exv::exv_method = exv;
            settings::form_factor::max_types = max_types;
            settings::form_factor::min_fraction = min_fraction;
            settings::hist::weighted_bins = weighted_bins;
        }
    };

    std::vector<double> q_values() {
        std::vector<double> q;
        for (double qi = 1e-3; qi < 5; qi *= 1.1) {q.push_back(qi);} // also far beyond the default q axis
        return q;
    }
}

TEST_CASE("form_factor_debye_transform: agrees with the direct sum") {
    settings::general::verbose = false;
    SettingsGuard guard;
    auto radiation = GENERATE(Radiation::XRay, Radiation::Neutron);
    settings::scattering::radiation = radiation;

    auto s = golomb_scatterers();
    auto molecule = as_molecule(s);
    auto q = q_values();
    auto expected = direct_sum(s, q);
    auto norm = scale(s, q);

    SECTION("default settings") {
        auto I = hist::form_factor_debye_transform(molecule, q);
        REQUIRE(I.size() == q.size());
        for (unsigned int k = 0; k < q.size(); ++k) {
            CHECK_THAT(I[k], Catch::Matchers::WithinAbs(expected[k], 1e-5*norm[k]));
        }
    }

    // the selection limits and the excluded volume model would otherwise fold the rare and the exotic types onto OTHER
    SECTION("restrictive form factor settings are lifted and restored") {
        settings::form_factor::max_types = 4;
        settings::form_factor::min_fraction = 0.5;
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;

        auto I = hist::form_factor_debye_transform(molecule, q);
        for (unsigned int k = 0; k < q.size(); ++k) {
            CHECK_THAT(I[k], Catch::Matchers::WithinAbs(expected[k], 1e-5*norm[k]));
        }
        CHECK(settings::form_factor::max_types == 4);
        CHECK(settings::form_factor::min_fraction == 0.5);
        CHECK(settings::exv::exv_method == settings::exv::ExvMethod::Fraser);
    }

    SECTION("the atomic weights are not used") {
        auto I = hist::form_factor_debye_transform(molecule, q);
        for (auto& body : molecule.get_bodies()) {
            for (auto& atom : body.get_atoms()) {atom.weight() = 1000;}
        }
        auto I_weighted = hist::form_factor_debye_transform(molecule, q);
        for (unsigned int k = 0; k < q.size(); ++k) {
            CHECK(I_weighted[k] == I[k]);
        }
    }
}

// on the default q axis, the calculation is the atomic and hydration part of the excluded volume models
TEST_CASE("form_factor_debye_transform: agrees with the average model without excluded volume") {
    settings::general::verbose = false;
    SettingsGuard guard;
    settings::molecule::implicit_hydrogens = true;
    settings::hist::weighted_bins = true;
    settings::form_factor::max_types = form_factor::total_ff_count;
    settings::form_factor::min_fraction = 0;
    settings::exv::exv_method = settings::exv::ExvMethod::Average;

    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();
    auto h = molecule.get_histogram();
    auto* model = dynamic_cast<hist::IScatteringModelExv*>(h.get());
    REQUIRE(model != nullptr);
    model->apply_excluded_volume_scaling_factor(0);
    auto expected = model->debye_transform();

    auto q = constants::axes::q_axis.sub_axis_covering(settings::axes::qmin, settings::axes::qmax).as_vector();
    auto I = hist::form_factor_debye_transform(molecule, q);
    REQUIRE(static_cast<int>(I.size()) == expected.size());
    for (unsigned int k = 0; k < q.size(); ++k) {
        CHECK_THAT(I[k], Catch::Matchers::WithinRel(expected[k], 1e-6));
    }
}

TEST_CASE("form_factor_debye_transform: atoms without form factors are refused") {
    settings::general::verbose = false;
    Molecule molecule({Body{std::vector<Atom>{Atom({0, 0, 0}, 1), Atom({1, 0, 0}, 1)}}});
    CHECK_THROWS_AS(hist::form_factor_debye_transform(molecule, {0.1, 0.2}), except::invalid_argument);
}
