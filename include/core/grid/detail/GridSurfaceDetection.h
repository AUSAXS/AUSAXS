// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <grid/detail/RadialLineGenerator.h>
#include <grid/detail/GridExcludedVolume.h>
#include <grid/GridFwd.h>

namespace ausaxs::grid::detail {
    class GridSurfaceDetection : private RadialLineGenerator {
        public:
            GridSurfaceDetection(observer_ptr<grid::Grid> grid);
            ~GridSurfaceDetection() override;

            exv::GridExcludedVolume detect() const;
            exv::GridExcludedVolume no_detect() const;

        private:
            observer_ptr<grid::Grid> grid;

            template<bool detect, bool unity_width>
            exv::GridExcludedVolume helper() const;

            bool collision_check(const Vector3<int>& loc) const;

            bool vacuum_collision_check(const Vector3<int>& loc) const;

            /**
             * @brief Determine the vacuum holes in the protein.
             *        This simply fills all non-water small gaps in the interior with vacuum voxels.
             */
            std::vector<Vector3<double>> determine_vacuum_holes() const;
    };
}