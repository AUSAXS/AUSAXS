#pragma once

#include <io/ExistingFile.h>

#include <string>
#include <vector>

namespace settings {
    namespace detail {
        void parse_option(const std::string& name, const std::vector<std::string>& value);
        bool is_comment_char(char c);
    }

    /**
     * @brief Read the settings from a file.
     */
    void read(const ::io::ExistingFile& path);

    /**
     * @brief Write the settings to a file.
     */
    void write(const ::io::File& path);

    /**
     * @brief Check if a settings file exists in the given directory, and read it if so.
     * @return True if a settings file was found and read.
     */
    bool discover(const ::io::Folder& path);
}