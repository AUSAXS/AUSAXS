// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/SimpleExvModel.h>
#include <hist/detail/CompactCoordinates.h>
#include <data/Molecule.h>
#include <utility/Logging.h>

using namespace ausaxs::hist::detail;

bool flag_simple_excluded_volume = false;
void SimpleExvModel::enable() {
    flag_simple_excluded_volume = true;
    logging::log("SimpleExvModel enabled.");
}

void SimpleExvModel::disable() {
    flag_simple_excluded_volume = false;
    logging::log("SimpleExvModel disabled.");
}

template<bool vbw>
void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates<vbw>& data_a, observer_ptr<const data::Molecule> molecule) {
    assert(molecule != nullptr && "SimpleExvModel::apply_simple_excluded_volume: molecule is nullptr.");
    if (flag_simple_excluded_volume) {
        assert(0 < molecule->size_atom() && "SimpleExvModel::apply_simple_excluded_volume: Division by zero. The molecule has no atoms.");
        data_a.implicit_excluded_volume(molecule->get_volume_grid()/molecule->size_atom());
    }
}
template void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates<true>& data_a, observer_ptr<const data::Molecule> molecule);
template void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates<false>& data_a, observer_ptr<const data::Molecule> molecule);