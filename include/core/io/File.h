#pragma once

#include <io/Folder.h>
#include <utility/TypeTraits.h>

#include <string_view>
#include <iosfwd>

namespace io {
    class File {
        public:
            File() = default;
            File(const io::Folder& folder, std::string_view name, std::string_view extension);
            File(std::string_view path);
            File(const char* path);
            File(const std::string& path);

            [[nodiscard]] std::string path() const;

            [[nodiscard]] std::string absolute_path() const;

            [[nodiscard]] operator std::string() const;

            /**
             * @brief Replace the extension of the file.
             */
            File& replace_extension(std::string_view extension) noexcept;

            /**
             * @brief Append to the name of the file.
             */
            File& append(std::string_view name) noexcept;

            /**
             * @brief Append to the name of the file.
             */
            [[nodiscard]] File append(std::string_view name) const noexcept;

            /**
             * @brief Get the stem of the file.
             */
            [[nodiscard]] std::string stem() const noexcept;

            /**
             * @brief Get the directory of the file.
             */
            [[nodiscard]] Folder& directory() noexcept;

            /**
             * @brief Get the directory of the file.
             */
            [[nodiscard]] const Folder& directory() const noexcept;

            /**
             * @brief Get the extension of the file, including the dot.
             */
            [[nodiscard]] std::string& extension() noexcept;

            /**
             * @brief Get the extension of the file, including the dot. 
             */
            [[nodiscard]] const std::string& extension() const noexcept;

            /**
             * @brief Create this file with the given contents.
             */
            void create(std::string_view contents = "") const; 

            /**
             * @brief Remove this file.
             */
            void remove() const;

            /**
             * @brief Check if the file exists.
             */
            [[nodiscard]] bool exists() const noexcept;

            /**
             * @brief Split a path into a directory, filename and extension.
             */
            static std::tuple<std::string, std::string, std::string> split(std::string_view path);

        private:
            Folder dir;
            std::string name;
            std::string ext;
    };
    static_assert(supports_nothrow_move_v<File>, "File should support nothrow move semantics.");
}

std::string operator+(const std::string& str, const io::File& file);
std::string operator+(const io::File& file, const std::string& str);
std::string operator+(const char* str, const io::File& file);
std::string operator+(const io::File& file, const char* str);
std::ostream& operator<<(std::ostream& os, const io::File& file);