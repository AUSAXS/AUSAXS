#pragma once

#include <data/DataFwd.h>
#include <io/IOFwd.h>

#include <vector>

namespace rigidbody {
    struct BodySplitter {
        /**
         * @brief Constructor. 
         * 
         * @param input The path to the input data file. 
         */
        static data::Molecule split(const io::File& input, std::vector<int> splits);
    };
}