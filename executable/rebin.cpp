#include <utility/SimpleDataset.h>

#include <filesystem>
#include <iostream>

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

    // load the input measurement
    string mfile = argv[1];
    SimpleDataset data(mfile);
    data.rebin();
    data.save(mfile + ".rebin.txt");

    std::filesystem::rename(mfile, mfile + ".original.txt");
    std::filesystem::rename(mfile + ".rebin.txt", mfile);
}