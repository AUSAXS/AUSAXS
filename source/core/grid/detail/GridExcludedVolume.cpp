// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <grid/detail/GridExcludedVolume.h>
#include <data/atoms/AtomFF.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Logging.h>

using namespace ausaxs::grid::exv;
using namespace ausaxs::data;

bool GridExcludedVolume::has_surface() const {
    return !surface.empty();
}

void GridExcludedVolume::save(const io::File& file) const {
    if (interior.empty() && surface.empty()) {
        logging::log("GridExcludedVolume::save: No interior or surface atoms to save. Skipping write.");
        return;
    }
    std::vector<AtomFF> atoms1, atoms2;
    
    for (int i = 0; i < static_cast<int>(interior.size()); ++i) {
        atoms1.emplace_back(interior[i], form_factor::form_factor_t::C);
    }

    for (int j = 0; j < static_cast<int>(surface.size()); ++j) {
        atoms2.emplace_back(surface[j], form_factor::form_factor_t::C);
    }

    std::vector<Body> bodies = {Body{atoms1}, Body{atoms2}};
    data::Molecule(bodies).save(file);
}