#pragma once

#include <utility/files/Folder.h>

namespace gmx {
    namespace detail {
        struct File {
            File();
            File(const std::string& path);
            File(const std::string& path, const std::string& ext);

            void validate(const std::string& path);

            operator std::string() const;

            bool empty() const;

            std::string parent_path() const;

            /**
             * @brief Delete this file from disk.
             */
            void remove() const;

            /**
             * @brief Create this file on disk.
             */
            void create(const std::string& content = "") const;

            /**
             * @brief Check if this file exists on disk.
             */
            bool exists() const;

            /**
             * @brief Move this file to another folder.
             */
            std::string move(const Folder& folder);

            /**
             * @brief Copy this file to another folder.
             * 
             * Any existing file with the same name will be overwritten.
             */
            virtual std::string copy(const Folder& folder) const;

            File& rename(const std::string& newname);

            /**
             * @brief Get the filename of this file.
             */
            std::string filename() const;

            /**
             * @brief Get the relative path to another file with this as the base.
             */
            std::string relative_path(const File& other) const;

            /**
             * @brief Get the relative path to another file with this as the base.
             */
            std::string relative_path(const std::string& s) const;

            /**
             * @brief Get the absolute path to this file.
             */
            std::string absolute() const;

            /**
             * @brief Get the stem of this file.
             * 
             * The stem is the filename without the extension.
             */
            std::string stem() const;

            std::string path;
            std::string ext;
        };
    }
}