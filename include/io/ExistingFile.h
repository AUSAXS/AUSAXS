#pragma once

#include <io/File.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile() = default;
            ExistingFile(const File& file);
            ExistingFile(File&& file);
            ExistingFile(const std::string& path);
            virtual ~ExistingFile() = default;

            ExistingFile& operator=(const std::string& path);

        private:
            void validate() const;
    };
}