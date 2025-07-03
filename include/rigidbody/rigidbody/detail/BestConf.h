// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <grid/GridFwd.h>
#include <data/atoms/Water.h>

#include <memory>
#include <vector>

namespace ausaxs::rigidbody::detail {
    struct BestConf {
        BestConf();
        BestConf(std::shared_ptr<grid::Grid> grid, std::vector<data::Water> waters, double chi2) noexcept;
        ~BestConf();

        std::shared_ptr<grid::Grid> grid;
        std::vector<data::Water> waters;
        double chi2;	
    };
}