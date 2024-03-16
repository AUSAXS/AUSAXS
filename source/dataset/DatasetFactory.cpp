/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <dataset/DatasetFactory.h>

#include <dataset/detail/DatasetReader.h>
#include <dataset/detail/DATReader.h>
#include <dataset/detail/XVGReader.h>
#include <constants/Constants.h>

std::unique_ptr<Dataset> factory::DatasetFactory::construct(const io::ExistingFile& file, unsigned int expected_cols) {
    std::unique_ptr<detail::DatasetReader> constructor;
    auto ext = utility::to_lowercase(file.extension());
    if (detail::DATReader::extensions.contains(ext)) {
        constructor = std::make_unique<detail::DATReader>();
    } else if (detail::XVGReader::extensions.contains(ext)) {
        constructor = std::make_unique<detail::XVGReader>();
    } else {
        throw except::invalid_operation("factory::create: Unknown file extension \"" + ext + "\".");
    }
    return constructor->construct(file, expected_cols);
}