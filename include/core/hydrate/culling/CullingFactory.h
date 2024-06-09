#pragma once

#include <hydrate/culling/CullingStrategy.h>
#include <settings/MoleculeSettings.h>
#include <data/DataFwd.h>

#include <memory>

namespace hydrate {
    namespace factory {
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> molecule);
        std::unique_ptr<CullingStrategy> construct_culling_strategy(observer_ptr<data::Molecule> grid, const settings::hydrate::CullingStrategy& choice);
    }
}