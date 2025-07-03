// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody {
    class Step {
        public: 
            Step(unsigned int n);

        private: 
            unsigned int iterations;
    };
}