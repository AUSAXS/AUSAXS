#include <CLI/CLI.hpp>

#include <data/Protein.h>

int main(int argc, char const *argv[]) {
    std::string pdb;
    CLI::App app{"Solvent density"};
    app.add_option("input", pdb, "The structure file.")->required();
    app.add_option("--output,-o", setting::general::output, "Path to save the generated figures at.")->default_val("output/solvent_density/");
    CLI11_PARSE(app, argc, argv);
    setting::general::output += utility::stem(pdb) + "/";

    Protein protein(pdb);
    

    return 0;
}