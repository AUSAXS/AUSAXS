#include <crystal/io/CrystalReaderFactory.h>
#include <crystal/io/GridReader.h>
#include <crystal/io/PDBReader.h>
#include <crystal/io/UnitCellReader.h>
#include <utility/Constants.h>
#include <utility/Exceptions.h>

#include <memory>

namespace crystal::io {
    std::unique_ptr<CrystalReader> CrystalReaderFactory::create(const std::string& filename) {
        if (constants::filetypes::unit_cell.validate(filename)) {
            return std::make_unique<UnitCellReader>();
        } else if (constants::filetypes::grid.validate(filename)) {
            return std::make_unique<GridReader>();
        } else if (constants::filetypes::structure.validate(filename)) {
            return std::make_unique<PDBReader>();
        } else {
            throw except::io_error("crystal::io::CrystalReaderFactory::create: Unknown file extension for file \"" + filename + "\"");
        }
    }
}