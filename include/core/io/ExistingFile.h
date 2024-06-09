#pragma once

#include <io/File.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile();
            ExistingFile(const File& file);
            ExistingFile(File&& file);
            ExistingFile(const std::string& path);
            ExistingFile(const char* path);
            virtual ~ExistingFile();

            ExistingFile& operator=(const std::string& path);

        private:
            void validate() const;
    };
}