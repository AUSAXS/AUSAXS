#pragma once

#include <utility/Axis.h>

/**
 * @brief Constants used to define the default axes.
 */
namespace ausaxs::constants::axes {
    using d_type = double;
    constexpr Axis d_axis(0, 2000, 8000);
    constexpr Axis q_axis(1e-4, 1, 200);
    constexpr auto q_vals = q_axis.as_array<q_axis.bins>(); 
    constexpr auto d_vals = d_axis.as_array<d_axis.bins>();
    constexpr double d_inv_width = 1./d_axis.width();
}