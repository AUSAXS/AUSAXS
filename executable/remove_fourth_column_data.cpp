#include <CLI/CLI.hpp>

#include <dataset/SimpleDataset.h>
#include <io/File.h>

int main(int argc, char const *argv[]) {
    io::File s_mfile;
    CLI::App app{"Remove the fourth column from a dataset."};
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    SimpleDataset dataset(s_mfile);
    dataset.save(s_mfile.append("_stripped"));
    return 0;
}