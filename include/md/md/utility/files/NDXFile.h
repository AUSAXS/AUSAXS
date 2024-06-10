#pragma once

#include <md/utility/files/File.h>

#include <fstream>

namespace md {
    // Index file
    struct NDXFile : public detail::File {
        NDXFile() = default;
        NDXFile(const std::string& name) : File(name, "ndx") {}
        NDXFile(const char* name) : NDXFile(std::string(name)) {}
        ~NDXFile() override = default;

        void append_file(const NDXFile& other) {
            std::ifstream ifs(other);
            std::ofstream ofs(*this, std::ios::app);
            ofs << ifs.rdbuf();
        }

        /**
         * @brief Check if the index file contains a group
         */
        bool contains(const std::string& group) const {
            std::ifstream ifs(*this);
            std::string line;
            while (std::getline(ifs, line)) {
                if (line.empty()) {continue;}
                if (line[0] != '[') {continue;}
                auto end = line.find(']');
                if (end == std::string::npos) {throw std::runtime_error("Invalid index file \"" + path + "\": line \"" + line + "\" does not contain a closing bracket.");}
                if (line[1] != ' ' || line[end-1] != ' ') {throw std::runtime_error("Invalid index file\"" + path + "\": line \"" + line + "\" does not contain a space before and after the group name.");}
                line = line.substr(2, end-3);
                if (line == group) {
                    return true;
                }
            }
            return false;
        }
    };
}