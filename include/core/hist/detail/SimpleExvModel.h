// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>
#include <utility/observer_ptr.h>

namespace ausaxs::hist::detail {
    /**
     * @brief The effective charge excluded volume approximation of settings::exv::ExvMethod::Simple, used by the non-form factor histogram managers.
     *        The managers apply it only when that method is selected; settings::exv::ExvMethod::None uses the same managers without it.
     */
    class SimpleExvModel {
        public:
            /**
             * @brief The excluded volume per atom of @a molecule, as consumed by apply_simple_excluded_volume.
             *
             * The excluded volume is that of the grid, which includes all symmetric copies, so it is shared over the atoms of all copies as well.
             * The grid is built and expanded lazily on first use, which is not thread-safe. This must therefore be called before any pool task
             * is submitted, never from inside one; the tasks only receive the result.
             *
             * @param molecule The molecule to use for the excluded volume calculation.
             */
            static double get_volume_per_atom(observer_ptr<const data::Molecule> molecule);

            /**
             * @brief Account for the excluded volume in the data.
             *		  Note: this should not be done for models with explicit excluded volume terms.
             *
             * This is done by subtracting the average excluded volume charge from each atom.
             *
             * @param data_a The atomic data to apply the excluded volume transformation to.
             * @param volume_per_atom The excluded volume per atom, see get_volume_per_atom.
             */
            static void apply_simple_excluded_volume(hist::detail::CompactCoordinates& data_a, double volume_per_atom);
    };
}