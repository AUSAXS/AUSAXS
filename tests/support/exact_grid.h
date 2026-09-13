#pragma once

#include <data/Body.h>
#include <grid/Grid.h>
#include <settings/GridSettings.h>
#include <utility/Limit3D.h>

#include <vector>

namespace ausaxs::test {
    /**
     * @brief A Grid with exactly the given axes.
     *
     * Grid itself only exposes content-sized constructors, which add the margin the hydration shell needs. Tests of the need to pin the 
     * axes instead so they can assert on specific bin coordinates, which is what this provides.
     */
    class ExactGrid : public grid::Grid {
        public:
            using Grid::operator=;

            explicit ExactGrid(const Limit3D& axes) : Grid(std::vector<data::Body>{}) {
                setup(Axis3D(axes, settings::grid::cell_width));
            }

            ~ExactGrid() override = default;
    };
}
