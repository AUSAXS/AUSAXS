#pragma once

#include <io/Folder.h>

#include <string>
#include <ostream>

namespace io {
    class File {
        public:
            File() = default;

            File(const io::Folder& folder, std::string_view name, std::string_view extension);

            File(const char* path);

            File(const std::string& path);

            virtual ~File() = default;

            [[nodiscard]] std::string path() const;

            void operator=(const std::string& path);

            [[nodiscard]] operator std::string() const;

            /**
             * @brief Replace the extension of the file.
             */
            void replace_extension(const std::string& extension) noexcept;

            /**
             * @brief Append to the name of the file.
             */
            void append(const std::string& name) noexcept;

            /**
             * @brief Append to the name of the file.
             */
            [[nodiscard]] File append(const std::string& name) const noexcept;

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
             * @brief Get the extension of the file.
             */
            [[nodiscard]] std::string& extension() noexcept;

            /**
             * @brief Get the extension of the file.
             */
            [[nodiscard]] const std::string& extension() const noexcept;

            /**
             * @brief Create this empty file.
             */
            void create() const; 

            /**
             * @brief Check if the file exists.
             */
            [[nodiscard]] bool exists() const noexcept;

            /**
             * @brief Split a path into a directory, filename and extension.
             */
            static std::tuple<std::string, std::string, std::string> split(const std::string& path);

        private:
            Folder dir;
            std::string name;
            std::string ext;
    };
}

std::string operator+(const std::string& str, const io::File& file);
std::string operator+(const io::File& file, const std::string& str);
std::string operator+(const char* str, const io::File& file);
std::string operator+(const io::File& file, const char* str);
std::ostream& operator<<(std::ostream& os, const io::File& file);
std::istringstream &operator>>(std::istringstream& in, io::File& val);