#pragma once

#include <grid/detail/RadialLineGenerator.h>
#include <grid/detail/GridExcludedVolume.h>
#include <grid/GridFwd.h>

namespace grid::detail {
    class GridSurfaceDetection : private RadialLineGenerator {
        public:
            GridSurfaceDetection(observer_ptr<grid::Grid> grid);
            ~GridSurfaceDetection() override;

            GridExcludedVolume detect_atoms();
            GridExcludedVolume detect_voxels();

        private:
            observer_ptr<grid::Grid> grid;

            bool collision_check(const Vector3<int>& loc);
    };
}