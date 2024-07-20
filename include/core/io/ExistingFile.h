#pragma once

#include <io/File.h>
#include <utility/TypeTraits.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile();
            ExistingFile(const File& file);
            ExistingFile(File&& file);

            template<::detail::string_like T>
            ExistingFile(const T& path) : ExistingFile(std::string_view(path)) {}
            ExistingFile(std::string_view path);

            template<::detail::string_like T>
            ExistingFile& operator=(const T& path) {return *this = std::string_view(path);}
            ExistingFile& operator=(std::string_view path);

        private:
            void validate() const;
    };
    static_assert(supports_nothrow_move_v<ExistingFile>, "ExistingFile should support nothrow move semantics.");
}