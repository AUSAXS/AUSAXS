// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/detail/IValidatedFile.h>
#include <stdexcept>

namespace ausaxs::io {
    namespace detail {
        struct validate_existing_file {
            static void validate(observer_ptr<File> f) {
                if (!f->empty() && !f->exists()) {
                    throw std::runtime_error("ExistingFile::validate: File \"" + f->path() + "\" does not exist.");
                }
            }
        };
    }

    class ExistingFile : public detail::IValidatedFile<detail::validate_existing_file> {
        using IValidatedFile::IValidatedFile;
    };
    static_assert(supports_nothrow_move_v<ExistingFile>, "ExistingFile should support nothrow move semantics.");
}