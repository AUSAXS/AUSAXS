// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

//? This is a poor place to put this
//? Consider moving all exv code in this dir out of grid to a new separate 'exv' folder
namespace ausaxs::grid::exv {
    /**
     * @brief Get the excluded volume of the currently selected excluded volume model. 
     * 
     * @param m The molecule to get the volume for.
     * @param d The excluded volume scaling factor.
     *
     * @return The excluded volume in cubic angstroms.
     */
    double get_volume_exv(observer_ptr<const data::Molecule> m, double d);
}