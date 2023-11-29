#pragma once

#include <utility/Limit.h>

/**
 * @brief A representation of a 2-dimensional limited span of values. 
 */
class Limit2D {
    public:
        Limit2D() noexcept;

        /**
         * @brief Construct a new 2D limit from two 1D limits. 
         * 
         * @param x The limit on the x-axis. 
         * @param y The limit on the y-axis. 
         */
        Limit2D(const Limit& x, const Limit& y) noexcept;

        /**
         * @brief Construct a new 2D limit from a minimum and maximum value along each axis. 
         */
        Limit2D(double xmin, double xmax, double ymin, double ymax) noexcept;

        /**
         * @brief Check if this object is fully initialized. Returns false if any of its Limits are empty.
         */
        [[nodiscard]] bool empty() const noexcept;

        Limit x; // The 1D limit on the x-axis. 
        Limit y; // The 1D limit on the y-axis.
};