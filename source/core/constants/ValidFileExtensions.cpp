// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <constants/ValidFileExtensions.h>
#include <utility/StringUtils.h>
#include <io/File.h>

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

std::string constants::filetypes::detail::guess_type(const io::File& path) {
    if (structure.check(path)) {
        return "structure";
    } if (saxs_data.check(path)) {
        return "saxs_data";
    } if (em_map.check(path)) {
        return "em_map";
    } if (unit_cell.check(path)) {
        return "unit_cell";
    } if (grid.check(path)) {
        return "grid";
    } if (config.check(path)) {
        return "config";
    }
    return "";
}

template class constants::filetypes::detail::FileType<1>;
template class constants::filetypes::detail::FileType<2>;
template class constants::filetypes::detail::FileType<3>;
template class constants::filetypes::detail::FileType<4>;