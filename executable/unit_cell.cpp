#include <utility/SimpleDataset.h>
#include <data/Protein.h>

#include <filesystem>
#include <iostream>

using std::string;

/**
 * @brief Create a CRYST1 record containing the entire structure. The original file is renamed. 
 */
int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;

    // check that we have at least one argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << "<pdbfile>" << std::endl;
        return 1;
    }

    // load the input file
    string pdbfile = argv[1];
    Protein protein(pdbfile);
    protein.generate_unit_cell();

    // save the new file
    // std::filesystem::rename(pdbfile, pdbfile + ".original.pdb");
    protein.save(pdbfile + "new.pdb");
}