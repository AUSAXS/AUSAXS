// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <constants/ValidFileExtensions.h>

#include <io/File.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

template<std::size_t N>
bool constants::filetypes::detail::FileType<N>::check(const io::File& path) const {
    if (!path.exists()) {return false;}
    auto file_ext = utility::to_lowercase(path.extension()); 
    for (const auto& ext : extensions) {
        if (file_ext == ext) {
            return true;
        }
    }
    return false;
}

template class constants::filetypes::detail::FileType<1>;
template class constants::filetypes::detail::FileType<2>;
template class constants::filetypes::detail::FileType<3>;
template class constants::filetypes::detail::FileType<4>;