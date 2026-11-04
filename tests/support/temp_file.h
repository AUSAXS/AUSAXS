// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/File.h>
#include <io/Folder.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace ausaxs::test {
    namespace detail {
        // A tag no other test file will carry. std::random_device is permitted to be deterministic, and is on some MinGW builds, so mix in the clock and a
        // counter: two test processes starting back-to-back must not agree on the names they hand out.
        inline std::string unique_tag() {
            static std::atomic<std::uint64_t> counter{0};
            std::uint64_t tag = std::random_device{}()
                ^ static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count())
                ^ (++counter << 48);

            std::ostringstream ss;
            ss << std::hex << tag;
            return ss.str();
        }
    }

    /**
     * @brief A scratch file under temp/, written on construction and deleted again when it goes out of scope.
     *
     * catch_discover_tests registers every TEST_CASE as its own CTest test, and CI runs them with `ctest --parallel`. Several processes of the *same* test
     * binary are therefore alive at once, each with its own copy of any function-local counter. A name built from a fixed prefix and a counter collides
     * across them: one process truncates and rewrites the file while another is parsing it, so a test silently runs someone else's input. The random tag
     * makes the name unique regardless of how many processes are in flight.
     */
    class TempFile : public io::File {
        public:
            /**
             * @brief Create a temporary file with the given extension and contents. 
             * @param extension File extension, including the dot.
             * @param contents The contents to write.
             */
            explicit TempFile(std::string_view extension, std::string_view contents = "")
                : io::File(io::Folder("temp"), "ausaxs_test_" + detail::unique_tag(), extension)
            {
                create(contents);
            }

            // move-only: a copy would share one path between two owners, and whichever is destroyed first deletes the file out from under the other.
            TempFile(const TempFile&) = delete;
            TempFile& operator=(const TempFile&) = delete;

            // NOLINTNEXTLINE(bugprone-use-after-move) - owns is a member of TempFile, not of the io::File base the mem-initialiser above moves from
            TempFile(TempFile&& other) noexcept : io::File(std::move(other)), owns(std::exchange(other.owns, false)) {}
            TempFile& operator=(TempFile&& other) noexcept {
                if (this != &other) {
                    if (owns) {remove();}
                    bool other_owns = std::exchange(other.owns, false);
                    io::File::operator=(std::move(other));
                    owns = other_owns;
                }
                return *this;
            }

            ~TempFile() override {if (owns) {remove();}}

        private:
            bool owns = true;
    };

    /**
     * @brief A scratch directory under temp/, created on construction and deleted with everything inside it when it goes out of scope.
     *
     * A CLI tool is handed an output folder and decides itself what to put there, so a test driving one cannot name the files it must clean up. The
     * random tag is there for the same reason as TempFile's: two parallel CTest processes of one binary must not share an output directory.
     */
    class TempFolder : public io::Folder {
        public:
            TempFolder() : io::Folder("temp/ausaxs_test_" + detail::unique_tag()) {create();}

            // move-only, as TempFile: a copy would share one path between two owners.
            TempFolder(const TempFolder&) = delete;
            TempFolder& operator=(const TempFolder&) = delete;

            // NOLINTNEXTLINE(bugprone-use-after-move) - owns is a member of TempFolder, not of the io::Folder base the mem-initialiser above moves from
            TempFolder(TempFolder&& other) noexcept : io::Folder(std::move(other)), owns(std::exchange(other.owns, false)) {}
            TempFolder& operator=(TempFolder&& other) noexcept {
                if (this != &other) {
                    if (owns) {discard();}
                    bool other_owns = std::exchange(other.owns, false);
                    io::Folder::operator=(std::move(other));
                    owns = other_owns;
                }
                return *this;
            }

            ~TempFolder() {if (owns) {discard();}}

        private:
            bool owns = true;

            // a destructor must not throw, and a folder the tool never wrote to is not an error
            void discard() const noexcept {
                std::error_code ec;
                std::filesystem::remove_all(path(), ec);
            }
    };
}
