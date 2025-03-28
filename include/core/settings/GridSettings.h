#pragma once

#include <utility/Limit3D.h>

namespace ausaxs::settings::grid {
    // The number of generated water molecules as a percent of the number of atoms.
    extern double water_scaling;

    // Cell width in Ångström. Each grid cell is a cube with sides of this length.
    extern double cell_width;

    // The percent increase in grid size in all dimensions when the grid size is automatically deduced based on an input vector of atoms.
    extern double scaling;
    
    // The radius of the excluded volume sphere around each atom. This does *not* block water placement and *only* affects volume calculations. 
    extern double min_exv_radius;

    // Whether to generate a cubic grid. This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations.
    extern bool   cubic;

    // The minimum number of bins in each dimension of the grid. 
    // This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations. A value of 0 means disabled. 
    extern unsigned int min_bins;

    // Grid-based excluded volume settings.
    namespace exv {
        // Whether to save the excluded volume grid when using the grid-based excluded volume calculations. This is primarily useful for debugging.
        extern bool   save;

        // The width of the excluded volume dummy atoms used for the grid-based excluded volume calculations in Å.
        extern double width;

        // The surface thickness of the grid in Ångström. This is used for fitting the excluded volume. 
        extern double surface_thickness;
    }

    namespace detail {
        extern double min_score; // min_score is the minimum percentage of radial lines which must not intersect anything to place a water molecule.
    }
}