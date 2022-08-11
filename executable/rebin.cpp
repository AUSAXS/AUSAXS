#include <utility/SimpleDataset.h>

#include <filesystem>
#include <iostream>
#include <fstream>

using std::string;

/**
 * @brief Rebin a dataset. This dramatically reduces the number of data points.
 */
int main(int argc, char const *argv[]) {
    // check that we have at least one argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << "<mfile>" << std::endl;
        return 1;
    }

    // check if file is already rebinned
    string mfile = argv[1];
    std::ifstream input(mfile);
    if (!input.is_open()) {
        std::cerr << "Error: Could not open file \"" << mfile << "\"" << std::endl;
        return 1;
    }
    string line;
    std::getline(input, line);
    if (line.find("REBINNED") != string::npos) {
        std::cerr << "Error: File \"" << mfile << "\" has already been rebinned" << std::endl;
        return 1;
    }

    // load the input measurement
    SimpleDataset data(mfile);
    data.rebin();
    data.save(mfile + ".rebin.txt", "REBINNED");

    std::filesystem::rename(mfile, mfile + ".original.txt");
    std::filesystem::rename(mfile + ".rebin.txt", mfile);
}