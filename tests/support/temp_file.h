// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/File.h>
#include <io/Folder.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include <string_view>

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
             * @param prefix File name prefix, without a directory. The unique tag is appended to it.
             * @param extension File extension, including the dot.
             * @param contents The contents to write.
             */
            TempFile(std::string_view prefix, std::string_view extension, std::string_view contents = "")
                : io::File(io::Folder("temp"), std::string(prefix) + "_" + detail::unique_tag(), extension)
            {
                create(contents);
            }

            ~TempFile() override {remove();}
    };
}
