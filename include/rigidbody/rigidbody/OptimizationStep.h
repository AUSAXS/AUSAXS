#pragma once

namespace ausaxs::rigidbody {
    class Step {
        public: 
            Step(unsigned int n);

        private: 
            unsigned int iterations;
    };
}