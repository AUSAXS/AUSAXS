#include <io/Writer.h>
#include <io/detail/PDBWriter.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

void io::Writer::write(const io::pdb::PDBStructure& s, const io::File& path) {
    auto ext = utility::to_lowercase(path.extension());
    if (ext == ".xml") { // .xml PDBStructure
        throw except::invalid_argument("PDBStructure::construct_writer: .xml input PDBStructures are not supported.");
    } else if (ext == ".pdb") { // .pdb PDBStructure
        io::detail::pdb::PDBWriter::write(s, path);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("PDBStructure::construct_writer: Unsupported extension \"" + path.extension() + "\" of input file \"" + path.str() + "\".");
    }
}