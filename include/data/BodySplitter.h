#pragma once

#include <algorithm>

#include <data/Protein.h>
#include <rigidbody/constraints/Constraint.h>
#include <io/ExistingFile.h>

namespace rigidbody {
    struct BodySplitter {
        /**
         * @brief Constructor. 
         * 
         * @param input The path to the input data file. 
         */
        static Protein split(const io::ExistingFile& input, std::vector<int> splits);
    };
}