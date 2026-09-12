// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hydrate/culling/CullingStrategy.h>
#include <settings/MoleculeSettings.h>

#include <memory>

namespace ausaxs::hydrate::factory {
    std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> molecule, bool global);
    std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> molecule, settings::hydrate::CullingStrategy choice);
}