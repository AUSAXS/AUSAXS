#pragma once

#include <io/Folder.h>
#include <string>

namespace io {
    class File {
        public:
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
             * @brief Get the stem of the file.
             */
            [[nodiscard]] std::string stem(const std::string& path);

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

std::string operator+(const char* str, const io::File& file);