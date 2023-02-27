#include <dataset/SimpleDataset.h>
#include <utility/Utility.h>

#include <filesystem>
#include <iostream>
#include <fstream>

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
    std::string mfile = argv[1];
    std::ifstream input(mfile);
    if (!input.is_open()) {
        std::cerr << "Error: Could not open file \"" << mfile << "\"" << std::endl;
        return 1;
    }
    std::string line;
    std::getline(input, line);
    if (line.find("REBINNED") != std::string::npos) {
        std::cerr << "Error: File \"" << mfile << "\" has already been rebinned" << std::endl;
        return 1;
    }

    // load the input measurement
    SimpleDataset data(mfile);
    data.rebin();

    auto path = std::filesystem::path(mfile).parent_path().string() + "/";
    data.save(path + "temp.dat", "REBINNED (prevents multiple rebinnings)");


    std::cout << path + utility::stem(mfile) + "_original.dat" << std::endl;
    std::filesystem::rename(mfile, path + utility::stem(mfile) + "_original.dat");
    std::filesystem::rename(path + "temp.dat", mfile);
}