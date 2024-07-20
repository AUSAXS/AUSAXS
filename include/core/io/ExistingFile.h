#pragma once

#include <io/detail/IValidatedFile.h>
#include <stdexcept>

namespace io {
    namespace detail {
        struct validate_existing_file {
            static void validate(observer_ptr<File> f) {
                if (!f->exists()) {
                    throw std::runtime_error("ExistingFile::validate: File \"" + f->path() + "\" does not exist.");
                }
            }
        };
    }

    class ExistingFile : public IValidatedFile<detail::validate_existing_file> {
        using IValidatedFile::IValidatedFile;
    };
    static_assert(supports_nothrow_move_v<ExistingFile>, "ExistingFile should support nothrow move semantics.");
}