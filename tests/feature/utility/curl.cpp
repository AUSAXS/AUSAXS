#include <catch2/catch_test_macros.hpp>

#include <utility/Curl.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>

#include <support/temp_file.h>

#include <fstream>

using namespace ausaxs;

TEST_CASE("Curl::download") {
    settings::general::verbose = false;
    test::TempFile file(".cif");
    curl::download("https://files.rcsb.org/ligands/view/LYS.cif", file);
    CHECK(file.exists());

    std::ifstream ifs(file);
    std::string line;
    std::getline(ifs, line);
    CHECK(line == "data_LYS");
}