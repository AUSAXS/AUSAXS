#include <crystal/io/CrystalReaderFactory.h>
#include <crystal/io/PDBReader.h>
#include <utility/Constants.h>
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