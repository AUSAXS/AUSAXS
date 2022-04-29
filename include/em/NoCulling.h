#pragma once

#include <em/CullingStrategy.h>

namespace em {
    class NoCulling : public CullingStrategy {
        public:
            NoCulling() {}

            ~NoCulling() override = default;

            std::vector<Atom> cull(std::list<Atom>& atoms) const override;
    };
}