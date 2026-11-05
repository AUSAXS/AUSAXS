// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/FormFactorDebyeCalculator.h>

#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <math/ConstexprMath.h>
#include <settings/ExvSettings.h>
#include <settings/FormFactorSettings.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <numeric>

using namespace ausaxs;

namespace {
    /**
     * @brief Exposes the form factor-resolved distance distributions of a molecule. The managers otherwise only hand them on to an
     *        excluded volume model, which is of no use here; the average one would even build the grid of the molecule for its
     *        excluded volume charge.
     */
    struct DistributionExporter : hist::HistogramManagerMTFFAvg<true> {
        using HistogramManagerMTFFAvg<true>::HistogramManagerMTFFAvg;
        auto distributions() {return this->compute_distributions();}
    };

    /**
     * @brief Activates a form factor set with a slot for every form factor type present in a molecule, by lifting the limits of
     *        settings::form_factor. The excluded volume model is disabled meanwhile, since the Fraser-based ones would fold every
     *        type without a tabulated excluded volume onto OTHER.
     *        The settings are restored afterwards; the active set is rebuilt with the restored excluded volume model.
     */
    struct EveryFormFactorScope {
        settings::exv::ExvMethod old_exv = settings::exv::exv_method;
        int old_max_types = settings::form_factor::max_types;
        double old_min_fraction = settings::form_factor::min_fraction;

        EveryFormFactorScope() {
            settings::exv::exv_method = settings::exv::ExvMethod::None;
            settings::form_factor::max_types = form_factor::total_ff_count;
            settings::form_factor::min_fraction = 0;
        }

        ~EveryFormFactorScope() {
            settings::form_factor::max_types = old_max_types;
            settings::form_factor::min_fraction = old_min_fraction;
            settings::exv::exv_method = old_exv;
        }

        EveryFormFactorScope(const EveryFormFactorScope&) = delete;
        EveryFormFactorScope& operator=(const EveryFormFactorScope&) = delete;
    };
}

std::vector<double> hist::form_factor_debye_transform(const data::Molecule& molecule, const std::vector<double>& q) {
    if (molecule.symmetry().has_symmetries()) {
        throw except::invalid_argument("hist::form_factor_debye_transform: Molecules with symmetries are not supported.");
    }
    for (const auto& atom : molecule.iterate_atoms()) {
        if (atom.form_factor_type() == form_factor::form_factor_t::UNKNOWN) {
            throw except::invalid_argument(
                "hist::form_factor_debye_transform: The molecule has atoms without form factor information, as happens when it is built from "
                "bare coordinates and weights. Form factors require the elements of the atoms, e.g. by loading the molecule from a structure file."
            );
        }
    }

    EveryFormFactorScope scope;
    form_factor::manager::use_form_factors(molecule);
    auto d = DistributionExporter(&molecule).distributions();
    Distribution3D<Shape::Triangular> p_aa(std::move(d.p_aa)); // each pair of distinct types in both orders, but only once per unordered pair of types
    Distribution2D p_aw(std::move(d.p_aw));                    // each atom-water pair once
    Distribution1D p_ww(std::move(d.p_ww));
    std::vector<double> d_axis = d.p_tot.get_weighted_axis();
    std::vector<double> I(q.size(), 0);
    if (d_axis.empty()) {return I;} // nothing to scatter
    d_axis[0] = 0; // as in DistanceHistogram: the first bin primarily holds the self-correlations, whose form factor is replaced below
    const int bins = static_cast<int>(d_axis.size());

    const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
    const int n_ff = form_factor::get_active_count();
    constexpr int first = form_factor::start_index_for_explicit_exv(); // the excluded volume slot is unused
    constexpr int water = form_factor::water_bin;

    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), [&] (int start, int end) {
        std::vector<double> sinqd(bins), f(n_ff), s(n_ff);
        auto debye_sum = [&sinqd] (auto begin, auto end) {return std::transform_reduce(begin, end, sinqd.begin(), 0.0);};
        for (int k = start; k < end; ++k) {
            for (int i = 0; i < bins; ++i) {sinqd[i] = math::sinc(q[k]*d_axis[i]);}
            for (int ff = first; ff < n_ff; ++ff) {
                auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
                f[ff] = form_factor::manager::evaluate_amplitude_at(type, q[k]);
                s[ff] = form_factor::manager::evaluate_self_at(type, q[k]);
            }

            double sum = 0;
            for (int ff1 = first; ff1 < n_ff; ++ff1) {
                for (int ff2 = ff1; ff2 < n_ff; ++ff2) {
                    sum += f[ff1]*f[ff2]*debye_sum(p_aa.begin(ff1, ff2), p_aa.end(ff1, ff2));
                }
                sum += 2*f[ff1]*f[water]*debye_sum(p_aw.begin(ff1), p_aw.end(ff1));

                // the first bin of a diagonal partial also holds the self-correlation of each atom, which the sum above evaluated with the
                // cross-correlation form factor; it is swapped for the self-correlation one. sinc(q*0) = 1
                sum += d.self_aa[ff1]*(s[ff1] - f[ff1]*f[ff1]);
            }
            sum += f[water]*f[water]*debye_sum(p_ww.begin(), p_ww.end());
            sum += d.self_ww*(s[water] - f[water]*f[water]);
            I[k] = sum;
        }
    });
    pool->wait();
    return I;
}
