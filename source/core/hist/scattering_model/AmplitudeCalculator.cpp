// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/AmplitudeCalculator.h>
#include <data/Body.h>  // IWYU pragma: keep
#include <data/Molecule.h>
#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/histogram_manager/detail/SimpleExvHelpers.h>
#include <math/ConstexprMath.h>
#include <math/Vector3.h>
#include <settings/ExvSettings.h>
#include <settings/InternalState.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <cmath>
#include <functional>

using namespace ausaxs;
using ausaxs::form_factor::form_factor_t;

namespace {
    /**
     * @brief A set of scatterers sharing one effective form factor shape. Scatterer j contributes weight_j*form_factor(q) to the amplitude.
     */
    struct ScattererClass {
        std::vector<Vector3<double>> positions{};
        std::vector<double> weights{};
        std::function<double(double)> form_factor;
        std::vector<double> distances{}; // from the centre, filled in by average_amplitude
    };

    double gaussian(double q) {return std::exp(-0.5*q*q);}

    /**
     * @brief The scatterers of the models without form factors, Simple and None.
     *        Their intensity is exp(-q^2) sum_ij w_i w_j sinc(q r_ij), so each scatterer is a point of weight w with the amplitude exp(-q^2/2).
     */
    std::vector<ScattererClass> weighted_scatterers(const data::Molecule& molecule, double displaced_charge) {
        ScattererClass atoms{.form_factor = gaussian};
        for (const auto& a : molecule.iterate_atoms()) {
            atoms.positions.emplace_back(a.coordinates());
            atoms.weights.emplace_back(a.weight() - displaced_charge);
        }

        ScattererClass waters{.form_factor = gaussian};
        for (const auto& w : molecule.iterate_waters()) {
            waters.positions.emplace_back(w.coordinates());
            waters.weights.emplace_back(w.weight());
        }

        std::vector<ScattererClass> classes;
        classes.emplace_back(std::move(atoms));
        classes.emplace_back(std::move(waters));
        return classes;
    }

    /**
     * @brief The scatterers of the form factor models, Vacuum and Fraser.
     *        The atoms are grouped by their active form factor slot, and every scatterer counts once, as in the form factor managers.
     *        With an explicit excluded volume, each atom also displaces the excluded volume form factor of its type.
     *        The waters carry no excluded volume of their own; see ScatteringModelExvBase.
     */
    std::vector<ScattererClass> form_factor_scatterers(const data::Molecule& molecule, bool explicit_exv) {
        form_factor::manager::use_form_factors(molecule);
        const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
        auto mapping = form_factor::manager::get_active_mapping();
        int n_ff = form_factor::get_active_count();

        auto exv_set = form_factor::ExvTableManager::get_current_exv_form_factor_set();
        std::vector<ScattererClass> classes(n_ff);
        for (int slot = 0; slot < n_ff; ++slot) {
            auto type = static_cast<form_factor_t>(ff_indices[slot]);
            bool has_exv = explicit_exv && type != form_factor_t::WATER && exv_set.contains(type);
            if (has_exv) {
                auto exv = exv_set.get(type);
                classes[slot].form_factor = [type, exv] (double q) {return form_factor::manager::evaluate(type, q) - exv.evaluate(q);};
            } else {
                classes[slot].form_factor = [type] (double q) {return form_factor::manager::evaluate(type, q);};
            }
        }

        for (const auto& a : molecule.iterate_atoms()) {
            if (a.form_factor_type() == form_factor_t::UNKNOWN) {
                throw except::invalid_argument(
                    "hist::average_amplitude: Attempted to use an atom with UNKNOWN form factor type.\n"
                    "Form factor information is required for the selected excluded volume model."
                );
            }
            auto& c = classes[mapping[static_cast<int>(a.form_factor_type())]];
            c.positions.emplace_back(a.coordinates());
            c.weights.emplace_back(1);
        }

        auto& waters = classes[form_factor::water_bin];
        for (const auto& w : molecule.iterate_waters()) {
            waters.positions.emplace_back(w.coordinates());
            waters.weights.emplace_back(1);
        }
        return classes;
    }

    std::vector<ScattererClass> get_scatterers(const data::Molecule& molecule) {
        using ExvMethod = settings::exv::ExvMethod;
        switch (settings::exv::exv_method) {
            case ExvMethod::None:
                return weighted_scatterers(molecule, 0);

            case ExvMethod::Simple: {
                if (molecule.size_atom() == 0) {return weighted_scatterers(molecule, 0);}
                double displaced_charge = settings::internal_state::solvent_density()*hist::detail::simple_exv::get_volume_per_atom(&molecule);
                return weighted_scatterers(molecule, displaced_charge);
            }

            case ExvMethod::Vacuum:
                return form_factor_scatterers(molecule, false);

            case ExvMethod::Fraser:
                return form_factor_scatterers(molecule, true);

            // we explicitly write each case to ensure we will get a compiler warning for new models in the future
            case ExvMethod::Average:
            case ExvMethod::Grid:
            case ExvMethod::GridSurface:
            case ExvMethod::GridScalable:
            case ExvMethod::CRYSOL:
            case ExvMethod::FoXS:
            case ExvMethod::Pepsi:
            case ExvMethod::WAXSiS:
                break;
        }
        throw except::invalid_argument(
            "hist::average_amplitude: The scattering amplitude is only available for the excluded volume models "
            "\"simple\", \"none\", \"vacuum\" and \"fraser\"."
        );
    }

    /**
     * @brief The centroid of the excess scattering length.
     *        If the excess nearly cancels, as for a contrast-matched particle, it is ill-defined, and the absolute excess is used instead.
     */
    Vector3<double> scattering_centre(const std::vector<ScattererClass>& classes) {
        Vector3<double> sum{0, 0, 0}, abs_sum{0, 0, 0};
        double total = 0, abs_total = 0;
        for (const auto& c : classes) {
            if (c.positions.empty()) {continue;}
            double f0 = c.form_factor(0);
            for (std::size_t i = 0; i < c.positions.size(); ++i) {
                double w = f0*c.weights[i];
                sum += c.positions[i]*w;
                abs_sum += c.positions[i]*std::abs(w);
                total += w;
                abs_total += std::abs(w);
            }
        }
        if (abs_total == 0) {throw except::invalid_argument("hist::average_amplitude: The molecule has no scattering length.");}
        if (std::abs(total) < 1e-6*abs_total) {return abs_sum/abs_total;}
        return sum/total;
    }
}

std::vector<double> hist::average_amplitude(const data::Molecule& molecule, const std::vector<double>& q_vals) {
    if (molecule.symmetry().has_symmetries()) {
        throw except::invalid_argument("hist::average_amplitude: Molecules with symmetries are not supported.");
    }

    auto classes = get_scatterers(molecule);
    auto centre = scattering_centre(classes);
    for (auto& c : classes) {
        c.distances.resize(c.positions.size());
        for (std::size_t i = 0; i < c.positions.size(); ++i) {c.distances[i] = c.positions[i].distance(centre);}
    }

    std::vector<double> A(q_vals.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q_vals.size()), [&classes, &q_vals, &A] (int start, int end) {
        for (int qi = start; qi < end; ++qi) {
            double q = q_vals[qi];
            for (const auto& c : classes) {
                if (c.positions.empty()) {continue;}
                double sum = 0;
                for (std::size_t i = 0; i < c.distances.size(); ++i) {sum += c.weights[i]*math::sinc(q*c.distances[i]);}
                A[qi] += c.form_factor(q)*sum;
            }
        }
    });
    pool->wait();
    return A;
}
