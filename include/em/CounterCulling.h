#pragma once

#include <em/CullingStrategy.h>

#include <vector>

namespace em {
    class CounterCulling : public CullingStrategy {
        public:
            CounterCulling() {}

            ~CounterCulling() override = default;

            std::vector<Atom> cull(std::list<Atom>& atoms) const override;
    };
}