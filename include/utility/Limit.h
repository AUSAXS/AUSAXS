#pragma once

#include <vector>
#include <array>
#include <initializer_list>
#include <math.h>
#include <ostream>
#include <string>

/**
 * @brief Limit
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
         * @brief Expand the limit by the specified percentage. 
         *        Both ends of the limit will be grow by the specified percentage. 
         */
        void expand(double percent) noexcept {
            double span = this->span();
            min -= span*percent;
            max += span*percent;
        } 


        Limit& operator-=(double c) noexcept {
            min-=c;
            max-=c;
            return *this;
        }
        Limit operator-(double c) const noexcept {return Limit(*this)-=c;}

        Limit& operator+=(double c) noexcept {
            min+=c;
            max+=c;
            return *this;
        }
        Limit operator+(double c) const noexcept {return Limit(*this)+=c;}

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