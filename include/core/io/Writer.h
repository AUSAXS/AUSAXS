#pragma once

#include <io/IOFwd.h>

namespace io::detail {
    /**
     * @brief Virtual super-class for all data file writers. 
     */
    class Writer {
        public:
            Writer() = default;
            Writer(const Writer&) = default;
            Writer(Writer&&) noexcept = default;
            Writer &operator=(const Writer&) = default;
            Writer &operator=(Writer&&) noexcept = default;
            virtual ~Writer() = default;

            bool operator==(const Writer& rhs) const = default;

            /**
             * @brief Write the contents of the backing File to a given path. 
             */
            virtual void write(const io::File&) = 0;
    };
}
