#pragma once

#include <utility/Limit3D.h>

namespace settings {
    namespace grid {
        extern double water_scaling; // The number of generated water molecules as a percent of the number of atoms.
        extern double width;         // Cell width in Ångström. Each grid cell is a cube with sides of this length.
        extern double scaling;       // The percent increase in grid size in all dimensions when the grid size is automatically deduced based on an input vector of atoms.
        extern double rvol;          // The radius of the excluded volume sphere around each atom. This does *not* block water placement and *only* affects volume calculations. 
        extern double exv_radius;    // The radius of the excluded volume sphere used for the grid-based excluded volume calculations in Å.
        extern bool   save_exv;      // Whether to save the excluded volume grid when using the grid-based excluded volume calculations. This is primarily useful for debugging.
        extern bool   cubic;         // Whether to generate a cubic grid. This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations.
        
        // The minimum number of bins in each dimension of the grid. 
        // This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations. A value of 0 means disabled. 
        extern unsigned int min_bins;

        namespace detail {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule.
        }
    }
}