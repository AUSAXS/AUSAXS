#pragma once

#include <utility/TypeTraits.h>
#include <io/IOFwd.h>

#include <string>
#include <vector>

namespace ausaxs::io {
    class Folder {
        public:
            Folder() = default;
            Folder(const Folder&) = default;
            Folder(Folder&&) noexcept = default;
            Folder& operator=(const Folder&) = default;
            Folder& operator=(Folder&&) noexcept = default;

            template<::detail::string_like T>
            Folder(const T& path) : Folder(std::string_view(path)) {}
            Folder(std::string_view path);

            /**
             * @brief Get the non-'/' terminated path of this folder.
             */
            [[nodiscard]] std::string path() const;

            /**
             * @brief Check if this folder exists on disk.
             */
            [[nodiscard]] bool exists() const noexcept;

            template<::detail::string_like T>
            void operator=(const T& path) {*this = std::string_view(path);}
            void operator=(std::string_view path);
            bool operator==(const Folder& other) const noexcept = default;

            [[nodiscard]] operator std::string() const;

            /**
             * @brief Get a list of all files inside this directory.
             */
            [[nodiscard]] std::vector<io::File> files() const;

            /**
             * @brief Get a list of all nested directories inside this directory.
             */
            [[nodiscard]] std::vector<io::Folder> directories() const;

            /**
             * @brief Check if this object is initialized.
             */
            [[nodiscard]] bool empty() const noexcept;

            /**
             * @brief Create this directory and all its parents.
             */
            void create() const;

        private:
            std::string dir;
    };
    static_assert(supports_nothrow_move_v<Folder>, "Folder should support nothrow move semantics.");
}

// Must return a string since the result can also be a file path 
[[nodiscard]] std::string operator+(std::string_view str, const io::Folder& folder);
[[nodiscard]] std::string operator+(const io::Folder& folder, std::string_view str);