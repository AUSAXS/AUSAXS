#include <catch2/catch_test_macros.hpp>

#include <settings/GeneralSettings.h>
#include <utility/Curl.h>

#include <support/temp_file.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

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

// the cases below need no network: a refused local port fails at once, and file:// URLs are served by libcurl itself
namespace {
    std::string read_all(const std::string& path) {
        std::ifstream ifs(path);
        return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
    }

    std::string file_url(const io::File& file) {
        return "file://" + std::filesystem::absolute(file.path()).generic_string();
    }
}

TEST_CASE("Curl::download: a failed transfer leaves no file behind") {
    settings::general::verbose = false;
    test::TempFile file(".cif");
    file.remove();

    SECTION("connection refused") {
        CHECK_FALSE(curl::download("http://127.0.0.1:1/LYS.cif", file));
    }
    SECTION("missing source") {
        test::TempFile missing(".cif");
        std::string url = file_url(missing);
        missing.remove();
        CHECK_FALSE(curl::download(url, file));
    }
    CHECK_FALSE(file.exists());
    CHECK_FALSE(std::filesystem::exists(file.path() + ".part"));
}

TEST_CASE("Curl::download: a failed transfer keeps an existing file") {
    settings::general::verbose = false;
    test::TempFile file(".cif", "data_LYS\n");
    CHECK_FALSE(curl::download("http://127.0.0.1:1/LYS.cif", file));
    CHECK(read_all(file.path()) == "data_LYS\n");
    CHECK_FALSE(std::filesystem::exists(file.path() + ".part"));
}

TEST_CASE("Curl::download: a successful transfer replaces the destination") {
    settings::general::verbose = false;
    test::TempFile source(".cif", "data_LYS\nnew contents\n");
    test::TempFile file(".cif", "stale");
    CHECK(curl::download(file_url(source), file));
    CHECK(read_all(file.path()) == "data_LYS\nnew contents\n");
    CHECK_FALSE(std::filesystem::exists(file.path() + ".part"));
}
