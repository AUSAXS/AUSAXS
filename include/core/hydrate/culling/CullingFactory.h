// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <settings/MoleculeSettings.h>
#include <data/DataFwd.h>

#include <memory>

namespace ausaxs::hydrate {
    namespace factory {
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> molecule, bool global);
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> molecule, settings::hydrate::CullingStrategy choice);
    }
}