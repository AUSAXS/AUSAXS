#pragma once

#include <io/IOFwd.h>

namespace io::detail {
    /**
     * @brief Virtual super-class for all data file readers. 
     */
    class Reader {
        public:
            Reader() = default;
            Reader(const Reader&) = default;
            Reader(Reader&&) noexcept = default;
            Reader &operator=(const Reader&) = default;
            Reader &operator=(Reader&&) noexcept = default;
            virtual ~Reader() = default;

            bool operator==(const Reader& rhs) const = default;

            /**
             * @brief Read the data stored in a file. 
             */
            virtual void read(const io::File&) = 0;
    };
}