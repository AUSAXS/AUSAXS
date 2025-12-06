// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <grid/exv/ExvVolume.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <utility/observer_ptr.h>
#include <settings/GridSettings.h>
#include <settings/ExvSettings.h>

double ausaxs::grid::exv::get_volume_exv(observer_ptr<const data::Molecule> m, double d) {
    assert(0 < m->size_atom() && "Molecule::get_volume_exv: Cannot compute excluded volume for an empty molecule.");
    auto fraser_helper = [m] () {
        double volume = 0;

        // we extract the volumes from the form factors since they have a better interface than the raw volume sets
        auto ff_table = form_factor::detail::ExvFormFactorSet(constants::exv::get_exv_set());
        for (const auto& body : m->get_bodies()) {
            volume += std::accumulate(body.get_atoms().begin(), body.get_atoms().end(), 0.0, [&ff_table] (double sum, const data::AtomFF& atom) {
                return sum + ff_table.get(atom.form_factor_type()).evaluate(0);
            });
        }
        return volume / constants::charge::density::water;
    };

    auto grid = m->get_grid();
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::Simple:
        case settings::exv::ExvMethod::Average:
            return m->get_volume_grid();
        
        case settings::exv::ExvMethod::Fraser:
            return fraser_helper()*hist::CompositeDistanceHistogramFFExplicit::exv_factor(0, d);

        case settings::exv::ExvMethod::Grid: 
        case settings::exv::ExvMethod::WAXSiS: {
            // note: not equivalent to grid volume! 
            // the grid can be finer than the resolution of the excluded volume, in which case every Nth bin is used
            auto exv_atoms = exv::RawGridExv::create(grid).interior.size();
            double single_vol = std::pow(settings::grid::cell_width, 3);
            return exv_atoms*single_vol;
        }

        case settings::exv::ExvMethod::GridScalable: {
            // scale the volume by the cubed factor
            auto exv = exv::RawGridExv::create(grid).interior.size();
            return exv*std::pow(settings::grid::cell_width*d, 3);
        }

        case settings::exv::ExvMethod::GridSurface: {
            // scale surface volumes by the factor
            auto exv = exv::RawGridWithSurfaceExv::create(grid);
            unsigned int interior_atoms = exv.interior.size();
            unsigned int exterior_atoms = exv.surface.size();
            double interior_vol = std::pow(settings::grid::cell_width, 3);
            double exterior_vol = std::pow(settings::grid::cell_width, 3)*hist::CompositeDistanceHistogramFFGridSurface::exv_factor(0, d);
            return interior_atoms*interior_vol + exterior_atoms*exterior_vol;
        }

        case settings::exv::ExvMethod::CRYSOL: {
            assert(m->size_atom() != 0 && "ExvVolume::get_volume_exv: Division by zero. The molecule has no atoms.");
            auto V = fraser_helper();
            return V*hist::CompositeDistanceHistogramCrysol::exv_factor(0, d, V/m->size_atom());
        }

        case settings::exv::ExvMethod::FoXS: {
            return fraser_helper()*hist::CompositeDistanceHistogramFoXS::exv_factor(0, d);
        }

        case settings::exv::ExvMethod::Pepsi: {
            return fraser_helper()*hist::CompositeDistanceHistogramPepsi::exv_factor(0, d);
        }

        default:
            throw std::runtime_error("Molecule::get_volume_exv: Unknown excluded volume method. Did you forget to add it to the switch statement?");
    }
}