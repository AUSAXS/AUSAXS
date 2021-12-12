#pragma once

#include "File.h"
#include "PDB_file.h"

#include <utility>

class FileConstructor {
    FileConstructor() = delete;

    static std::unique_ptr<File> construct(const string& path) {
        if (path.find(".xml") != string::npos || path.find(".XML") != string::npos) { // .xml file
            print_err("Error in Protein::Protein: .xml input files are not supported.");
            exit(1);
        } else if (path.find(".pdb") != string::npos || path.find(".PDB") != string::npos) { // .pdb file
            return std::make_unique<PDB_file>(path);
        } else { // anything else - we cannot handle this
            print_err((format("Error in Protein::Protein: Invalid file extension of input file %1%.") % path).str());
            exit(1);
        }
    }
};