#pragma once

#include <io/ExistingFile.h>

#include <string>

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