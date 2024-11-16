#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_pdb_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".pdb") {throw std::runtime_error("PDBFile::validate: File \"" + f->path() + "\" is not a protein data bank file (.pdb).");}
            }
        };
    }

    // Protein data bank file
    struct PDBFile : public io::detail::IValidatedFile<detail::validate_pdb_file> {
        using IValidatedFile::IValidatedFile;
    };
}