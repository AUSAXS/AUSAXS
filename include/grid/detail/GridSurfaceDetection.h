#pragma once

#include <grid/detail/RadialLineGenerator.h>
#include <grid/detail/GridExcludedVolume.h>
#include <grid/GridFwd.h>

namespace grid::detail {
    class GridSurfaceDetection : private RadialLineGenerator {
        public:
            GridSurfaceDetection(observer_ptr<grid::Grid> grid);
            ~GridSurfaceDetection() override;

            GridExcludedVolume detect() const;
            GridExcludedVolume no_detect() const;

        private:
            observer_ptr<grid::Grid> grid;

            template<bool detect>
            GridExcludedVolume helper() const;

            bool collision_check(const Vector3<int>& loc) const;
    };
}