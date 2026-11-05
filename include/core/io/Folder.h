// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <utility/TypeTraits.h>

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

            template<ausaxs::detail::string_like T>
            Folder(const T& path) : Folder(std::string_view(path)) {}
            Folder(std::string_view path);

            /**
             * @brief Get the path of this folder, without trailing separators.
             */
            [[nodiscard]] std::string path() const;
            [[nodiscard]] std::string str() const; //< @copydoc path

            /**
             * @brief Check if this folder exists on disk.
             */
            [[nodiscard]] bool exists() const noexcept;

            template<ausaxs::detail::string_like T>
            Folder& operator=(const T& path) {*this = std::string_view(path); return *this;}
            Folder& operator=(std::string_view path);
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

            /**
             * @brief Delete this directory and everything inside it. Does nothing if it does not exist.
             *
             * @throws std::filesystem::filesystem_error if it could not be removed.
             */
            void remove() const;

            /**
             * @brief Move this directory on disk, such that it becomes @a destination.
             *        A file is moved into a folder, but a folder is itself a location, so @a destination is where it ends up rather than what it
             *        is placed inside. @a destination must not already exist.
             *
             * @return The moved directory.
             * @throws std::filesystem::filesystem_error if it could not be moved.
             */
            io::Folder move(const io::Folder& destination) const;

        private:
            std::string dir;
    };
    static_assert(supports_nothrow_move_v<Folder>, "Folder should support nothrow move semantics.");
}