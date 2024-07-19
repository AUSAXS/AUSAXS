#pragma once

#include <io/File.h>
#include <utility/TypeTraits.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile();
            ExistingFile(const File& file);
            ExistingFile(File&& file);

            template<::detail::string_type T>
            ExistingFile(const T& path);

            template<::detail::string_type T>
            ExistingFile& operator=(const T& path);

        private:
            void validate() const;
    };
    static_assert(supports_nothrow_move_v<ExistingFile>, "ExistingFile should support nothrow move semantics.");
}