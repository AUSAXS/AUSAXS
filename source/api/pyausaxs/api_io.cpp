// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_io.h>
#include <constants/ValidFileExtensions.h>
#include <io/ExistingFile.h>

using namespace ausaxs;

API bool io_is_pdb(
    const char* path,
    int* status
) {return execute_with_catch([&]() {
    return constants::filetypes::structure.check(io::ExistingFile(path));
}, status);}

API bool io_is_saxs_data(
    const char* path,
    int* status
) {return execute_with_catch([&]() {
    return constants::filetypes::saxs_data.check(io::ExistingFile(path));
}, status);}

API bool io_is_em_map(
    const char* path,
    int* status
) {return execute_with_catch([&]() {
    return constants::filetypes::em_map.check(io::ExistingFile(path));
}, status);}

API bool io_is_rigidbody_config(
    const char* path,
    int* status
) {return execute_with_catch([&]() {
    return constants::filetypes::config.check(io::ExistingFile(path));
}, status);}