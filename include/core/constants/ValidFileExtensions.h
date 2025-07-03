// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>

#include <array>

namespace ausaxs::constants::filetypes {
    namespace detail {
        template<std::size_t N>
        class FileType {
            public:
                consteval FileType(std::array<const char*, N> extensions) : extensions(std::move(extensions)) {}

                /**
                * @brief Check if a file exists and has one of the allowed extensions for this file type.
                */
                bool check(const io::File& path) const;

            private:
                std::array<const char*, N> extensions;
        };
    }

    constexpr detail::FileType structure = {std::array{".pdb",    ".ent",  ".cif"         }};
    constexpr detail::FileType saxs_data = {std::array{".dat",    ".rsr",  ".xvg"         }};
    constexpr detail::FileType em_map    = {std::array{".map",    ".ccp4", ".mrc", ".rec" }};
    constexpr detail::FileType unit_cell = {std::array{".cell",   ".uc"                   }};
    constexpr detail::FileType grid      = {std::array{".grid"                            }};
    constexpr detail::FileType config    = {std::array{".config", ".conf", ".txt"         }};
    // * remember to correct cpp file if more than 4 possible extensions are added
}