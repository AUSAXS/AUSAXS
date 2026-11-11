// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/ManagerResults.h>

#include <data/Molecule.h>
#include <form_factor/FormFactorInfo.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/detail/BinEstimate.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/histogram_manager/detail/SimpleExvHelpers.h>
#include <hist/scattering_model/ScatteringModelAverage.h>
#include <hist/scattering_model/ScatteringModelFraser.h>
#include <hist/scattering_model/ScatteringModelNoExv.h>
#include <hist/scattering_model/ScatteringModelSimple.h>
#include <hist/scattering_model/crysol/ScatteringModelCrysol.h>
#include <hist/scattering_model/foxs/ScatteringModelFoXS.h>
#include <hist/scattering_model/pepsi/ScatteringModelPepsi.h>
#include <settings/ExvSettings.h>
#include <settings/InternalState.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

namespace {
    /**
     * @brief Construct a form factor scattering model from the exported distributions, followed by any model-specific arguments.
     */
    template<typename Model, bool wb, typename... Args>
    std::unique_ptr<IScatteringModel> make_ff_model(ManagerDistributions<wb, true>&& d, Args&&... args) { // NOLINT
        return std::make_unique<Model>(
            Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
            Distribution2D(std::move(d.p_aw)),
            Distribution1D(std::move(d.p_ww)),
            std::move(d.self_aa),
            d.self_ww,
            std::move(d.p_tot),
            std::forward<Args>(args)...
        );
    }
}

template<bool wb, bool ff>
ManagerDistributions<wb, ff> hist::detail::export_distributions(distance_calculator::HistogramStore<wb>& store, int aa, int aw, int ww) {
    ManagerDistributions<wb, ff> res;
    if constexpr (ff) { // read before the results are exported
        res.self_aa = store.get_self_3d(aa);
        res.self_ww = store.get_self_1d(ww);
    }
    res.p_ww = store.export_1d(ww);
    res.p_tot = typename ManagerDistributions<wb, ff>::ww_t(store.bins());
    auto add = [&p_tot = res.p_tot] (auto begin) {std::transform(p_tot.begin(), p_tot.end(), begin, p_tot.begin(), std::plus<>());};
    if constexpr (ff) {
        res.p_aa = store.export_3d(aa);
        res.p_aw = store.export_2d(aw);

        // no atom has the excluded volume type, so its rows are empty; they are skipped to make that explicit
        int n_ff = form_factor::get_active_count();
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_ff; ++ff1) {
            for (int ff2 = ff1; ff2 < n_ff; ++ff2) {add(res.p_aa.begin(ff1, ff2));}
            add(res.p_aw.begin(ff1));
        }
    } else {
        res.p_aa = store.export_1d(aa);
        res.p_aw = store.export_1d(aw);
        add(res.p_aa.begin());
        add(res.p_aw.begin());
    }
    add(res.p_ww.begin());

    // downsize our axes to only the relevant area
    res.resize(hist::detail::trimmed_bin_count(res.p_tot));
    return res;
}

template<bool wb>
std::unique_ptr<IScatteringModel> hist::detail::make_explicit_model(
    ManagerDistributions<wb, true>&& d, observer_ptr<const data::Molecule> protein // NOLINT
) {
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::Vacuum:  return make_ff_model<ScatteringModelNoExv> (std::move(d));
        case settings::exv::ExvMethod::FoXS:    return make_ff_model<ScatteringModelFoXS>  (std::move(d));
        case settings::exv::ExvMethod::Pepsi:   return make_ff_model<ScatteringModelPepsi> (std::move(d), protein);
        case settings::exv::ExvMethod::CRYSOL:  return make_ff_model<ScatteringModelCrysol>(std::move(d), protein);
        default:                                return make_ff_model<ScatteringModelFraser>(std::move(d));
    }
}

template<bool wb>
std::unique_ptr<IScatteringModel> hist::detail::make_average_model(
    ManagerDistributions<wb, true>&& d, observer_ptr<const data::Molecule> protein // NOLINT
) {
    double Z_exv_avg = protein->size_atom() == 0 ? 0 : protein->get_volume_grid()*settings::internal_state::solvent_density()/protein->size_atom();
    return make_ff_model<ScatteringModelAverage>(std::move(d), Z_exv_avg);
}

template<bool wb>
std::unique_ptr<IScatteringModel> hist::detail::make_simple_model(
    ManagerDistributions<wb, true>&& d, observer_ptr<const data::Molecule> protein // NOLINT
) {
    // the plain managers weigh each atom by its charge less the displaced charge per atom. All atoms of a type share both, so the same
    // weights can be applied to the histogram of each type instead
    double displaced_charge = settings::internal_state::solvent_density()*simple_exv::get_volume_per_atom(protein);
    const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
    int n_ff = form_factor::get_active_count();
    std::vector<double> w(n_ff);
    for (int ff = form_factor::start_index_for_explicit_exv(); ff < n_ff; ++ff) {
        w[ff] = constants::charge::get_ff_charge(static_cast<form_factor::form_factor_t>(ff_indices[ff])) - displaced_charge;
    }
    double w_water = constants::charge::get_ff_charge(static_cast<form_factor::form_factor_t>(ff_indices[form_factor::water_bin]));

    Distribution3D<hist::Shape::Triangular> p_aa_ff(std::move(d.p_aa));
    Distribution2D p_aw_ff(std::move(d.p_aw));
    int bins = static_cast<int>(d.p_tot.size());
    Distribution1D p_aa(bins), p_aw(bins), p_ww(std::move(d.p_ww));
    for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_ff; ++ff1) {
        for (int ff2 = ff1; ff2 < n_ff; ++ff2) {
            for (int i = 0; i < bins; ++i) {p_aa.index(i) += w[ff1]*w[ff2]*p_aa_ff.index(ff1, ff2, i);}
        }

        // the form factor histograms count each atom-water pair once, the plain ones in both orders
        for (int i = 0; i < bins; ++i) {p_aw.index(i) += 2*w[ff1]*w_water*p_aw_ff.index(ff1, i);}
    }
    for (int i = 0; i < bins; ++i) {p_ww.index(i) *= w_water*w_water;}

    // the total of the form factor managers is a count; it only lends its axis and bin centres, and the values are rebuilt from the parts
    auto h = std::make_unique<ScatteringModelSimple>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(d.p_tot));
    h->apply_water_scaling_factor(1);
    return h;
}

template<bool wb, bool ff>
std::unique_ptr<IScatteringModel> hist::detail::make_model(
    ManagerDistributions<wb, ff>&& d, observer_ptr<const data::Molecule> protein // NOLINT
) {
    if constexpr (!ff) {
        return std::make_unique<ScatteringModelSimple>(
            Distribution1D(std::move(d.p_aa)),
            Distribution1D(std::move(d.p_aw)),
            Distribution1D(std::move(d.p_ww)),
            std::move(d.p_tot)
        );
    } else {
        if (settings::exv::exv_method == settings::exv::ExvMethod::Simple) {return make_simple_model<wb>(std::move(d), protein);}
        if (settings::exv::exv_method == settings::exv::ExvMethod::Average) {return make_average_model<wb>(std::move(d), protein);}
        return make_explicit_model<wb>(std::move(d), protein);
    }
}

template ManagerDistributions<false, false> hist::detail::export_distributions<false, false>(distance_calculator::HistogramStore<false>&, int, int, int);
template ManagerDistributions<false, true> hist::detail::export_distributions<false, true>(distance_calculator::HistogramStore<false>&, int, int, int);
template ManagerDistributions<true, false> hist::detail::export_distributions<true, false>(distance_calculator::HistogramStore<true>&, int, int, int);
template ManagerDistributions<true, true> hist::detail::export_distributions<true, true>(distance_calculator::HistogramStore<true>&, int, int, int);
template std::unique_ptr<IScatteringModel> hist::detail::make_model<false, false>(ManagerDistributions<false, false>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_model<false, true>(ManagerDistributions<false, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_model<true, false>(ManagerDistributions<true, false>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_model<true, true>(ManagerDistributions<true, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_explicit_model<false>(ManagerDistributions<false, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_explicit_model<true>(ManagerDistributions<true, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_average_model<false>(ManagerDistributions<false, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_average_model<true>(ManagerDistributions<true, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_simple_model<false>(ManagerDistributions<false, true>&&, observer_ptr<const data::Molecule>);
template std::unique_ptr<IScatteringModel> hist::detail::make_simple_model<true>(ManagerDistributions<true, true>&&, observer_ptr<const data::Molecule>);
