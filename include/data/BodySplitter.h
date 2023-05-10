#pragma once

#include <vector>

class Protein;
namespace io {class ExistingFile;}
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