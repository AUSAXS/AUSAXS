// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/DatasetFactory.h>

#include <constants/ValidFileExtensions.h>
#include <io/detail/dataset/DATReader.h>
#include <io/detail/dataset/DatasetReader.h>
#include <io/detail/dataset/XVGReader.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

std::unique_ptr<Dataset> factory::DatasetFactory::construct(const io::ExistingFile& file, int expected_cols) {
    std::unique_ptr<io::detail::dataset::DatasetReader> constructor;
    auto ext = utility::to_lowercase(file.extension());
    if (io::detail::dataset::DATReader::extensions.contains(ext)) {
        constructor = std::make_unique<io::detail::dataset::DATReader>();
    } else if (io::detail::dataset::XVGReader::extensions.contains(ext)) {
        constructor = std::make_unique<io::detail::dataset::XVGReader>();
    } else {
        if (auto format = constants::filetypes::detail::guess_type(file); !format.empty()) {
            throw except::invalid_operation("DatasetFactory::construct: Expected a data file format, but got what appears to be a " + format + " file.");
        }
        throw except::invalid_operation("factory::create: Unknown file extension \"" + ext + "\".");
    }
    return constructor->construct(file, expected_cols);
}