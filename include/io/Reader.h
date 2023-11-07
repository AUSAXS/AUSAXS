#pragma once

#include <io/IOFwd.h>

#include <string>

namespace io::detail {
    /**
     * @brief Virtual super-class for all data file readers. 
     */
    class Reader {
        public:
            virtual ~Reader() = default;

            /**
             * @brief Read the data stored in a file. 
             */
            virtual void read(const io::File&) = 0;

            bool operator==(const Reader& rhs) const = default;
    };
}