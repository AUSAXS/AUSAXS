#pragma once

#include <io/File.h>
#include <utility/TypeTraits.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile();
            ExistingFile(const File& file);
            ExistingFile(File&& file);
            ExistingFile(const std::string& path);
            ExistingFile(const char* path);

            ExistingFile& operator=(const std::string& path);

        private:
            void validate() const;
    };
    static_assert(supports_nothrow_move_v<ExistingFile>, "ExistingFile should support nothrow move semantics.");
}