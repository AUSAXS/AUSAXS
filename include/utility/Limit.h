#pragma once

#include <vector>
#include <array>
#include <initializer_list>
#include <math.h>

/**
 * @brief \class Limit
 * 
 * A simple representation of a limited span of values. 
 */
class Limit {
    public:
        /**
         * @brief Default constructor.
         */
        Limit() noexcept : min(0), max(0) {}

        /**
         * @brief Constructor. 
         * 
         * @param min Minimum value of the limit. 
         * @param max Maximum value of the limit. 
         */
        Limit(double min, double max) noexcept : min(min), max(max) {}

        /**
         * @brief Get the interval spanned by this limit. 
         */
        double span() const noexcept {return max-min;}

        /**
         * @brief Get the center of this limit.
         */
        double center() const noexcept {return min + span()/2;}

        /**
         * @brief Merge two limits, returning their combined span.
         */
        void merge(const Limit& other) noexcept {
            min = std::min(min, other.min);
            max = std::max(max, other.max);
        }

        /**
         * @brief Equality operator. Check if this Limit is equal to another.
         */
        bool operator==(const Limit& rhs) const noexcept {return min == rhs.min && max == rhs.max;}

        /**
         * @brief Inequality operator.
         * 
         * Check if this object is different from another. 
         */
        bool operator!=(const Limit& rhs) const noexcept {return !operator==(rhs);}

        /**
         * @brief Subtract a constant from both ends of this limit. 
         */
        Limit operator-(double rhs) const noexcept {return Limit(min-rhs, max-rhs);}

        /**
         * @brief Stream output operator. 
         * 
         * Allows this object to easily be output to a given stream. 
         */
        friend std::ostream& operator<<(std::ostream& os, const Limit& axis) noexcept {os << axis.to_string(); return os;}

        /**
         * @brief Get a string representation of this object. 
         */
        std::string to_string(std::string prepend = "") const noexcept {return prepend + "(" + std::to_string(min) + ", " + std::to_string(max) + ")";}

        /**
         * @brief Check if this Limit is empty (min == max == 0).
         */
        [[nodiscard]] 
        bool empty() const noexcept {return min == 0 && max == 0;}

        double min; // The minimum value of this limit. 
        double max; // The maximum value of this limit. 
};

/**
 * @brief \class Limit3D
 * 
 * A representation of a 3-dimensional limited span of values. 
 */
class Limit3D {
    public:
        Limit3D() noexcept {}
        /**
         * @brief Constructor.
         * 
         * Construct a new 3D limit based on three 1D limits. 
         * 
         * @param x The limit on the x-axis. 
         * @param y The limit on the y-axis. 
         * @param z The limit on the z-axis. 
         */
        Limit3D(const Limit& x, const Limit& y, const Limit& z) noexcept : x(x), y(y), z(z) {}

        /**
         * @brief Constructor. 
         * 
         * Construct a new 3D limit based on the minimum and maximum values along each axis. 
         */
        Limit3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) noexcept : x(xmin, xmax), y(ymin, ymax), z(zmin, zmax) {}

        /**
         * @brief Check if this object is fully initialized. Returns false if any of its Limits are empty.
         */
        [[nodiscard]] 
        bool empty() const noexcept {return x.empty() || y.empty() || z.empty();}

        Limit x; // The 1D limit on the x-axis. 
        Limit y; // The 1D limit on the y-axis. 
        Limit z; // The 1D limit on the z-axis. 
};