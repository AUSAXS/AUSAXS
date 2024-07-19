#pragma once

#include <io/Folder.h>
#include <utility/TypeTraits.h>

#include <string_view>
#include <iosfwd>

namespace io {
    class File {
        public:
            File() = default;
            File(const File&) = default;
            File(File&&) noexcept = default;
            File &operator=(const File&) = default;
            File &operator=(File&&) noexcept = default;
            File(const io::Folder& folder, std::string_view name, std::string_view extension);
            File(std::string_view name, std::string_view extension);

            template<::detail::string_type T>
            File(const T& path);

            virtual ~File() = default;

            /**
             * @brief Get the path to this file relative to the current directory.
             */
            [[nodiscard]] std::string path() const;

            /**
             * @brief Get the absolute path to this file.
             */
            [[nodiscard]] std::string absolute_path() const;

            /**
             * @brief Get the relative path to another file with this as the base.
             */
            [[nodiscard]] std::string relative_path(const File& other) const;

            [[nodiscard]] operator std::string() const;

            /**
             * @brief Replace the extension of the file.
             */
            File& replace_extension(std::string_view extension) noexcept;

            /**
             * @brief Append to the name of the file.
             */
            File& append(std::string_view name) noexcept;
            [[nodiscard]] File append(std::string_view name) const noexcept; //< @copydoc append()

            /**
             * @brief Get the filename.
             *       This is just the combined stem and extension of the file.
             */
            [[nodiscard]] std::string filename() const noexcept;

            /**
             * @brief Get the stem of the file.
             */
            [[nodiscard]] std::string& stem() noexcept;
            [[nodiscard]] std::string stem() const noexcept; //< @copydoc stem()

            /**
             * @brief Get the directory of the file.
             */
            [[nodiscard]] Folder& directory() noexcept;
            [[nodiscard]] const Folder& directory() const noexcept; //< @copydoc directory()

            /**
             * @brief Get the extension of the file, including the dot.
             */
            [[nodiscard]] std::string& extension() noexcept;
            [[nodiscard]] const std::string& extension() const noexcept; //< @copydoc extension()

            /**
             * @brief Create this file on disk with the given contents.
             *        Parent directories are created if they do not already exist.
             *        If the file already exists, it is overwritten.
             */
            void create(std::string_view contents = "") const; 

            /**
             * @brief Delete this file from disk. 
             */
            void remove() const;

            /**
             * @brief Move this file to the given folder.
             *        The destination folder is created if it does not already exist.
             *
             * @returns The path to the new file.
             */
            io::File move(const io::Folder& folder) const;

            /**
             * @brief Copy this file to the given folder.
             *        The destination folder is created if it does not already exist.
             *
             * @returns The path to the new file.
             */
            io::File copy(const io::Folder& folder) const;

            /**
             * @brief Rename this file on disk.
             *
             * @returns The path to the new file.
             */
            io::File rename(std::string_view name) const;

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