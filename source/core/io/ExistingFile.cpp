/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/ExistingFile.h>
#include <utility/Exceptions.h>

io::ExistingFile::ExistingFile() = default;

io::ExistingFile::ExistingFile(const File& file) : File(file) {
    validate();
}

io::ExistingFile::ExistingFile(File&& file) : File(std::move(file)) {
    validate();
}

template<::detail::string_type T>
io::ExistingFile::ExistingFile(const T& path) : File(path) {
    validate();
}

template<::detail::string_type T>
io::ExistingFile& io::ExistingFile::operator=(const T& path) {
    File::operator=(path);
    validate();
    return *this;
}

void io::ExistingFile::validate() const {
    if (!exists()) {
        throw except::io_error("io::ExistingFile: File \"" + path() + "\" does not exist.");
    }
}

template io::ExistingFile::ExistingFile(const char* const&);
template io::ExistingFile::ExistingFile(const std::string&);
template io::ExistingFile::ExistingFile(const std::string_view&);
template io::ExistingFile& io::ExistingFile::operator=(const char* const&);
template io::ExistingFile& io::ExistingFile::operator=(const std::string&);
template io::ExistingFile& io::ExistingFile::operator=(const std::string_view&);