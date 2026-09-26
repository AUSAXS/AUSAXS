// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/ManagerResults.h>

#include <constants/Constants.h>
#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool wb, bool ff>
ManagerDistributions<wb, ff> hist::detail::export_distributions(distance_calculator::HistogramStore<wb>& store, int aa, int aw, int ww) {
    ManagerDistributions<wb, ff> res;
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
std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_explicit_histogram(
    ManagerDistributions<wb, true>&& d, settings::exv::ExvMethod method, observer_ptr<const data::Molecule> protein
) {
    switch (method) {
        case settings::exv::ExvMethod::FoXS:
            return std::make_unique<CompositeDistanceHistogramFoXS>(
                Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
                Distribution2D(std::move(d.p_aw)),
                Distribution1D(std::move(d.p_ww)),
                std::move(d.p_tot)
            );
        case settings::exv::ExvMethod::Pepsi:
            return std::make_unique<CompositeDistanceHistogramPepsi>(
                Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
                Distribution2D(std::move(d.p_aw)),
                Distribution1D(std::move(d.p_ww)),
                std::move(d.p_tot),
                protein
            );
        case settings::exv::ExvMethod::CRYSOL:
            return std::make_unique<CompositeDistanceHistogramCrysol>(
                Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
                Distribution2D(std::move(d.p_aw)),
                Distribution1D(std::move(d.p_ww)),
                std::move(d.p_tot),
                protein
            );
        default:
            return std::make_unique<CompositeDistanceHistogramFFExplicit>(
                Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
                Distribution2D(std::move(d.p_aw)),
                Distribution1D(std::move(d.p_ww)),
                std::move(d.p_tot)
            );
    }
}

template<bool wb, bool ff>
std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_histogram(
    ManagerDistributions<wb, ff>&& d, settings::exv::ExvMethod method, observer_ptr<const data::Molecule> protein
) {
    if constexpr (!ff) {
        return std::make_unique<CompositeDistanceHistogram>(
            Distribution1D(std::move(d.p_aa)),
            Distribution1D(std::move(d.p_aw)),
            Distribution1D(std::move(d.p_ww)),
            std::move(d.p_tot)
        );
    } else {
        if (method != settings::exv::ExvMethod::Average) {return make_explicit_histogram<wb>(std::move(d), method, protein);}
        double Z_exv_avg = protein->size_atom() == 0 ? 0 : protein->get_volume_grid()*constants::charge::density::water/protein->size_atom();
        return std::make_unique<CompositeDistanceHistogramFFAvg>(
            Distribution3D<hist::Shape::Triangular>(std::move(d.p_aa)),
            Distribution2D(std::move(d.p_aw)),
            Distribution1D(std::move(d.p_ww)),
            std::move(d.p_tot),
            Z_exv_avg
        );
    }
}

template ManagerDistributions<false, false> hist::detail::export_distributions<false, false>(distance_calculator::HistogramStore<false>&, int, int, int);
template ManagerDistributions<false, true> hist::detail::export_distributions<false, true>(distance_calculator::HistogramStore<false>&, int, int, int);
template ManagerDistributions<true, false> hist::detail::export_distributions<true, false>(distance_calculator::HistogramStore<true>&, int, int, int);
template ManagerDistributions<true, true> hist::detail::export_distributions<true, true>(distance_calculator::HistogramStore<true>&, int, int, int);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_histogram<false, false>(ManagerDistributions<false, false>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_histogram<false, true>(ManagerDistributions<false, true>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_histogram<true, false>(ManagerDistributions<true, false>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_histogram<true, true>(ManagerDistributions<true, true>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_explicit_histogram<false>(ManagerDistributions<false, true>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
template std::unique_ptr<ICompositeDistanceHistogram> hist::detail::make_explicit_histogram<true>(ManagerDistributions<true, true>&&, settings::exv::ExvMethod, observer_ptr<const data::Molecule>);
