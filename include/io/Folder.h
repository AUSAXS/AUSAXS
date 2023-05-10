#pragma once

#include <string>
#include <vector>

namespace io {
    class Folder {
        public:
            Folder();

            Folder(const std::string& path);

            virtual ~Folder();

            [[nodiscard]] std::string path() const;

            [[nodiscard]] bool exists() const noexcept;

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
}

std::string operator+(const char* str, const io::Folder& folder);
std::string operator+(const io::Folder& folder, const char* str);
std::string operator+(const std::string& str, const io::Folder& folder);
std::string operator+(const io::Folder& folder, const std::string& str);