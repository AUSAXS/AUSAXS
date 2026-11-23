// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/AmplitudeFormFactors.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/AmplitudeHistogram.h>
#include <hist/scattering_model/ScatteringModelFraser.h>
#include <settings/ExvSettings.h>
#include <utility/Exceptions.h>

#include <cmath>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    /**
     * @brief The form factors of ScatteringModelSimple: every scatterer is a point of its weight with the amplitude exp(-q^2/2), since the
     *        intensity is exp(-q^2) sum_ij w_i w_j sinc(q r_ij). The waters are scaled by c_w.
     */
    class GaussianFormFactors : public AmplitudeFormFactors {
        public:
            std::vector<double> evaluate(const ScatteringParameters& params, double q) const override {
                double g = std::exp(-0.5*q*q);
                std::vector<double> res(2);
                res[AmplitudeScatterers::atom_class] = g;
                res[AmplitudeScatterers::water_class] = params.cw*g;
                return res;
            }

            bool form_factors() const override {return false;}
    };

    /**
     * @brief The form factors of ScatteringModelNoExv: B_a(q)*f_t(q) for the atoms of each form factor slot, with the atomic Debye-Waller
     *        factor B_a, and c_w*f_w(q) for the waters.
     */
    class VacuumFormFactors : public AmplitudeFormFactors {
        public:
            std::vector<double> evaluate(const ScatteringParameters& params, double q) const override {
                const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
                const int n_ff = form_factor::get_active_count();
                const double B = ScatteringModelFraser::get_atomic_debye_waller_factor(q, params.DW_sigma_atomic);
                constexpr int water = form_factor::water_bin;
                std::vector<double> g(n_ff, 0); // no scatterer has the excluded volume slot
                for (int ff = form_factor::start_index_for_explicit_exv(); ff < n_ff; ++ff) {
                    if (ff == water) {continue;}
                    g[ff] = B*form_factor::manager::evaluate(static_cast<form_factor::form_factor_t>(ff_indices[ff]), q);
                }
                g[water] = params.cw*form_factor::manager::evaluate(static_cast<form_factor::form_factor_t>(ff_indices[water]), q);
                return g;
            }

            bool form_factors() const override {return true;}
    };

    /**
     * @brief The form factors of ScatteringModelFraser: those of VacuumFormFactors, less the excluded volume of each atomic type scaled as in
     *        the intensity, rho*cx(q)*B_x(q)*fx_t(q). The waters carry no excluded volume, see ScatteringModelExvBase.
     */
    class FraserFormFactors : public VacuumFormFactors {
        public:
            std::vector<double> evaluate(const ScatteringParameters& params, double q) const override {
                auto g = VacuumFormFactors::evaluate(params, q);

                // the excluded volume tables of the intensity are only tabulated on the default q axis, so its form factors are evaluated directly
                const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
                const auto exv_set = form_factor::ExvTableManager::get_current_exv_form_factor_set();
                const double x = params.crho*ScatteringModelFraser::exv_factor(q, params.cx)
                    *ScatteringModelFraser::get_exv_debye_waller_factor(q, params.DW_sigma_exv);
                for (int ff = form_factor::start_index_for_explicit_exv(); ff < form_factor::get_active_count(); ++ff) {
                    if (ff == form_factor::water_bin) {continue;}
                    auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
                    if (!exv_set.contains(type)) {continue;} // a type without a volume has no excluded volume in the intensity either
                    g[ff] -= x*exv_set.get(type).evaluate(q);
                }
                return g;
            }
    };
}

std::shared_ptr<const AmplitudeFormFactors> hist::make_amplitude_form_factors() {
    using ExvMethod = settings::exv::ExvMethod;
    switch (settings::exv::exv_method) {
        case ExvMethod::None:
        case ExvMethod::Simple:
            return std::make_shared<GaussianFormFactors>();
        case ExvMethod::Vacuum:
            return std::make_shared<VacuumFormFactors>();
        case ExvMethod::Fraser:
            return std::make_shared<FraserFormFactors>();

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
    return nullptr;
}

std::vector<double> hist::amplitude_transform(
    const AmplitudeFormFactors* form_factors, const AmplitudeScatterers& scatterers, const std::vector<double>& q,
    const ScatteringParameters& params, const std::optional<Vector3<double>>& centre
) {
    if (form_factors == nullptr) {
        throw except::invalid_argument(
            "hist::amplitude_transform: The scattering amplitude is only available for the excluded volume models "
            "\"simple\", \"none\", \"vacuum\" and \"fraser\"."
        );
    }
    if (form_factors->form_factors() != scatterers.form_factors) {
        throw except::invalid_argument(
            "hist::amplitude_transform: The scatterers do not have the layout of the model. "
            "Were they taken from a histogram manager built for another excluded volume model?"
        );
    }

    auto g = [form_factors, &params] (double q) {return form_factors->evaluate(params, q);};
    auto origin = centre.has_value() ? *centre : AmplitudeHistogram::scattering_centre(scatterers, g(0));
    return AmplitudeHistogram(scatterers, origin).transform(q, g);
}
