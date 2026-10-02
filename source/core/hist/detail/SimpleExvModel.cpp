// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/SimpleExvModel.h>

#include <data/Molecule.h>
#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <hist/detail/CompactCoordinates.h>

using namespace ausaxs::hist::detail;

void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates& data_a, observer_ptr<const data::Molecule> molecule) {
    assert(molecule != nullptr && "SimpleExvModel::apply_simple_excluded_volume: molecule is nullptr.");
    int n_atoms = molecule->symmetry().size_atom_total();
    assert(0 < n_atoms && "SimpleExvModel::apply_simple_excluded_volume: Division by zero. The molecule has no atoms.");
    data_a.implicit_excluded_volume(molecule->get_volume_grid()/static_cast<double>(n_atoms));
}