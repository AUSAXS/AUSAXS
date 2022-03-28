#pragma once

#include <em/CullingStrategy.h>

namespace em {
    class CounterCulling : public CullingStrategy {
        public:
            CounterCulling() {}

            ~CounterCulling() override = default;

            vector<Atom> cull(list<Atom>& atoms) const override;
    };
}