// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <vector>
#include <array>
#include <initializer_list>
#include <math.h>
#include <ostream>
#include <string>

namespace ausaxs {
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
            Limit() noexcept;

            /**
             * @brief Constructor. 
             * 
             * @param min Minimum value of the limit. 
             * @param max Maximum value of the limit. 
             */
            Limit(double min, double max) noexcept;

            /**
             * @brief Get the interval spanned by this limit. 
             */
            double span() const noexcept;

            /**
             * @brief Get the center of this limit.
             */
            double center() const noexcept;

            /**
             * @brief Merge two limits, returning their combined span.
             */
            void merge(const Limit& other) noexcept;

            /**
             * @brief Expand the limit by the specified percentage. 
             *        Both ends of the limit will be grow by the specified percentage. 
             */
            void expand(double percent) noexcept;


            Limit& operator-=(double c) noexcept;

            Limit operator-(double c) const noexcept;

            Limit& operator+=(double c) noexcept;

            Limit operator+(double c) const noexcept;

            /**
             * @brief Equality operator. Check if this Limit is equal to another.
             */
            bool operator==(const Limit& rhs) const noexcept;

            /**
             * @brief Inequality operator.
             * 
             * Check if this object is different from another. 
             */
            bool operator!=(const Limit& rhs) const noexcept;

            /**
             * @brief Stream output operator. 
             * 
             * Allows this object to easily be output to a given stream. 
             */
            friend std::ostream& operator<<(std::ostream& os, const Limit& axis) noexcept {os << axis.to_string(); return os;}

            /**
             * @brief Get a string representation of this object. 
             */
            std::string to_string(const std::string& prepend = "") const noexcept;

            /**
             * @brief Check if this Limit is empty (min == max == 0).
             */
            [[nodiscard]] bool empty() const noexcept;

            double min; // The minimum value of this limit. 
            double max; // The maximum value of this limit. 
    };
}