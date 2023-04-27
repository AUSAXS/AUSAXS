#pragma once

#include <io/File.h>

#include <string>

/**
 * @brief Virtual super-class for all data file writers. 
 */
class Writer {
    public:
        virtual ~Writer() = default;

        /**
         * @brief Write the contents of the backing File to a given path. 
         */
        virtual void write(const io::File&) = 0;
};