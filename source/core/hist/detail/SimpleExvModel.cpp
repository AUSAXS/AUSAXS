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

void SimpleExvModel::apply_simple_excluded_volume(hist::detail::CompactCoordinates& data_a, observer_ptr<const data::Molecule> protein) {
    if (flag_simple_excluded_volume) {
        data_a.implicit_excluded_volume(protein->get_volume_grid()/protein->size_atom());
    }
}