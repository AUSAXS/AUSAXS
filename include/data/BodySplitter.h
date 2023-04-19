#pragma once

#include <algorithm>

#include <data/Protein.h>
#include <rigidbody/constraints/Constraint.h>

namespace rigidbody {
    struct BodySplitter {
        /**
         * @brief Constructor. 
         * 
         * @param input The path to the input data file. 
         */
        static Protein split(const std::string& input, std::vector<int> splits);
    };
}