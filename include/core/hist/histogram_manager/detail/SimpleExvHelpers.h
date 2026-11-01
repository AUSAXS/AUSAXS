// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>
#include <utility/observer_ptr.h>

/**
 * @brief The effective charge excluded volume approximation of settings::exv::ExvMethod::Simple.
 */
namespace ausaxs::hist::detail::simple_exv {
	/**
	 * @brief The excluded volume per atom of @a molecule, as consumed by apply.
	 *
	 * @param molecule The molecule to use for the excluded volume calculation.
	 */
	double get_volume_per_atom(observer_ptr<const data::Molecule> molecule);

	/**
	 * @brief Account for the excluded volume in the data.
	 *        This is done by subtracting the average excluded volume charge from each atom.
	 *
	 * @param data_a The atomic data to apply the excluded volume transformation to.
	 * @param volume_per_atom The excluded volume per atom, see get_volume_per_atom.
	 */
	void apply(hist::detail::CompactCoordinates& data_a, double volume_per_atom);
}
