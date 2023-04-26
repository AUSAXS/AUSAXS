#pragma once

#include <string>
#include <vector>

namespace io {
    class Folder {
        public:
            Folder() = default;

            Folder(const std::string& path);

            virtual ~Folder() = default;

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