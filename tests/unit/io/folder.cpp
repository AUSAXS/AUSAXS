#include <catch2/catch_test_macros.hpp>

#include <io/Folder.h>

#include <filesystem>
#include <fstream>

using namespace ausaxs;

TEST_CASE("Folder::Folder") {
    SECTION("empty") {
        io::Folder folder("");
        CHECK(folder.path() == ".");
    }

    SECTION("simple") {
        io::Folder folder("test");
        CHECK(folder.path() == "test");
    }

    SECTION("simple with trailing slash") {
        io::Folder folder("tests/");
        CHECK(folder.path() == "tests");
    }

    SECTION("nested") {
        io::Folder folder("tests/folder");
        CHECK(folder.path() == "tests/folder");
    }

    SECTION("nested with trailing slash") {
        io::Folder folder("tests/folder/");
        CHECK(folder.path() == "tests/folder");
    }

    SECTION("repeated trailing slashes") {
        io::Folder folder("tests/folder//");
        CHECK(folder.path() == "tests/folder");
    }

    SECTION("root") {
        io::Folder folder("/");
        CHECK(folder.path() == "/");
        CHECK_FALSE(folder.empty());
    }
}

TEST_CASE("Folder::exists") {
    SECTION("false") {
        io::Folder folder("fake");
        CHECK(folder.exists() == false);

        io::Folder folder2("fake/folder");
        CHECK(folder2.exists() == false);
    }

    SECTION("true") {
        io::Folder folder("tests");
        CHECK(folder.exists());

        io::Folder folder2("tests/files");
        CHECK(folder2.exists());
    }
}

TEST_CASE("Folder::create") {
    SECTION("simple") {
        std::string path = "___dummy";
        io::Folder folder(path);
        folder.create();
        CHECK(folder.exists());
        std::filesystem::remove(path);
    }

    SECTION("nested") {
        std::string path = "dummy/nested";
        io::Folder folder(path);
        folder.create();
        CHECK(folder.exists());
        std::filesystem::remove(path);
    }
}

TEST_CASE("Folder::empty") {
    SECTION("true") {
        io::Folder folder;
        CHECK(folder.empty());
    }

    SECTION("false") {
        io::Folder folder("test");
        CHECK_FALSE(folder.empty());
    }
}

TEST_CASE("Folder::str") {
    io::Folder folder("tests");
    CHECK(folder.str() == "tests");
}

TEST_CASE("Folder::operator==") {
    io::Folder f1("tests");
    io::Folder f2("tests");
    io::Folder f3("other");
    CHECK(f1 == f2);
    CHECK_FALSE(f1 == f3);
}

TEST_CASE("Folder::remove") {
    std::filesystem::path root = "temp/tests/folder_remove";
    std::filesystem::remove_all(root);

    SECTION("with contents") {
        io::Folder folder((root/"a").string());
        io::Folder((root/"a"/"b").string()).create();
        std::ofstream((root/"a"/"b"/"file.txt").string()) << "content";
        folder.remove();
        CHECK_FALSE(folder.exists());
        CHECK(io::Folder(root.string()).exists());
    }

    SECTION("missing") {
        io::Folder folder((root/"missing").string());
        CHECK_NOTHROW(folder.remove());
    }

    SECTION("uninitialized") {
        CHECK_NOTHROW(io::Folder().remove());
    }
    std::filesystem::remove_all(root);
}

TEST_CASE("Folder::move") {
    std::filesystem::path root = "temp/tests/folder_move";
    std::filesystem::remove_all(root);

    io::Folder from((root/"from").string());
    io::Folder((root/"from"/"nested").string()).create();
    io::Folder to((root/"to").string());

    auto moved = from.move(to);
    CHECK(moved == to);
    CHECK_FALSE(from.exists());
    CHECK(io::Folder((root/"to"/"nested").string()).exists());
    std::filesystem::remove_all(root);
}
