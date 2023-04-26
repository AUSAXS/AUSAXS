#pragma once

#include <io/File.h>

namespace io {
    class ExistingFile : public File {
        public:
            ExistingFile(const std::string& path);
            virtual ~ExistingFile() = default;
    };
}