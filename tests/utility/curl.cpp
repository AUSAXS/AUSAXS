#include <catch2/catch_test_macros.hpp>

#include <utility/Curl.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>

#include <fstream>

TEST_CASE("Curl::download") {
    settings::general::verbose = false;
    io::File file("temp/test/curl/LYS.cif");
    file.create();
    curl::download("https://files.rcsb.org/ligands/view/LYS.cif", file);
    CHECK(file.exists());

    std::ifstream ifs(file);
    std::string line;
    std::getline(ifs, line);
    CHECK(line == "data_LYS");
}