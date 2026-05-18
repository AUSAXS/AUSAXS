// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs {
    /**
     * @brief Abstract base class for matrix decompositions.
     *
     * Derived classes implement decompose() to factorise a matrix; the shared @c precision constant
     * defines the tolerance used when comparing values during the decomposition.
     */
    class Decomposition {
        public: 
            Decomposition() {}
            virtual ~Decomposition() {}
            virtual void decompose() = 0;

        protected:
            static constexpr double precision = 1e-9;
    };
}