#pragma once

#include <em/CullingStrategy.h>

namespace em {
    class NoCulling : public CullingStrategy {
        public:
            NoCulling() {}

            ~NoCulling() override = default;

            vector<Atom> cull(list<Atom>& atoms) const override;
    };
}