#include <io/ExistingFile.h>
#include <utility/Exceptions.h>

io::ExistingFile::ExistingFile(const std::string& path) : File(path) {
    if (!exists()) {
        throw except::io_error("io::ExistingFile: File \"" + path + "\" does not exist.");
    }
}
