#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_ndx_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".ndx") {throw std::runtime_error("NDXFile::validate: File \"" + f->path() + "\" is not an index file (.ndx).");}
            }
        };
    }

    // Index file
    struct NDXFile : public io::detail::IValidatedFile<detail::validate_ndx_file> {
        using IValidatedFile::IValidatedFile;

        void append_file(const NDXFile& other);

        /**
         * @brief Check if the index file contains a group
         */
        bool contains(const std::string& group) const;
    };
}