#pragma once

#include <string>

namespace io {class ExistingFile;}
/**
 * @brief Virtual super-class for all data file readers. 
 */
class Reader {
    public:
        virtual ~Reader() = default;

        /**
         * @brief Read the data stored in a file. 
         */
        virtual void read(const io::ExistingFile&) = 0;
};