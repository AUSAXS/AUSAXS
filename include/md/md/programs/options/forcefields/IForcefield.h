#pragma once

#include <string>
#include <memory>

namespace ausaxs::md::option {
    enum class Forcefield {
        AMBER99SB,
        AMBER99SB_ILDN,
    };

    struct IForcefield {
        static std::unique_ptr<IForcefield> construct(Forcefield ff);

        virtual std::string filename() const = 0;
        bool exists() const;
    };
}