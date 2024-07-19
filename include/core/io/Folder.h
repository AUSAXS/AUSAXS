#pragma once

#include <utility/TypeTraits.h>

#include <string>
#include <vector>

namespace io {
    class Folder {
        public:
            Folder() = default;
            Folder(const std::string& path);
            Folder(const char* path);

            [[nodiscard]] std::string path() const;

            [[nodiscard]] bool exists() const noexcept;

            [[nodiscard]] Folder& operator+(const std::string& str) noexcept;

            [[nodiscard]] Folder& operator+(const Folder& folder) noexcept;

            void operator=(const std::string& path);

            [[nodiscard]] operator std::string() const;

            [[nodiscard]] std::vector<std::string> files() const;

            [[nodiscard]] std::vector<std::string> directories() const;

            /**
             * @brief Create this directory and all its parents.
             */
            void create() const;

        private:
            std::string dir;
    };
    static_assert(supports_nothrow_move_v<Folder>, "Folder should support nothrow move semantics.");
}

std::string operator+(const char* str, const io::Folder& folder);
std::string operator+(const std::string& str, const io::Folder& folder);
