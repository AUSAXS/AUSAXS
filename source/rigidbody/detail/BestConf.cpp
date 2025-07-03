// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/detail/BestConf.h>

using namespace ausaxs;

rigidbody::detail::BestConf::BestConf() = default;
rigidbody::detail::BestConf::BestConf(std::shared_ptr<grid::Grid> grid, std::vector<data::Water> waters, double chi2) noexcept : grid(std::move(grid)), waters(std::move(waters)), chi2(chi2) {}
rigidbody::detail::BestConf::~BestConf() = default;

