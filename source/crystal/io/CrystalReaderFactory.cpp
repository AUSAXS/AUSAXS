/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <crystal/io/CrystalReaderFactory.h>
#include <crystal/io/PDBReader.h>
#include <crystal/io/CrystalReader.h>
#include <crystal/io/GridReader.h>
#include <crystal/io/UnitCellReader.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <utility/Exceptions.h>

#include <memory>

namespace crystal::factory {
    std::unique_ptr<crystal::io::CrystalReader> CrystalReaderFactory::create(const ::io::ExistingFile& filename) {
        if (constants::filetypes::unit_cell.validate(filename)) {
            return std::make_unique<crystal::io::UnitCellReader>();
        } else if (constants::filetypes::grid.validate(filename)) {
            return std::make_unique<crystal::io::GridReader>();
        } else if (constants::filetypes::structure.validate(filename)) {
            return std::make_unique<crystal::io::PDBReader>();
        } else {
            throw except::io_error("crystal::io::CrystalReaderFactory::create: Unknown file extension for file \"" + filename + "\"");
        }
    }
}