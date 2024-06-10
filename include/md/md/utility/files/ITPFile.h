#pragma once

#include <md/utility/files/File.h>

#include <vector>

namespace gmx {
    // Include topology file
    struct ITPFile : public detail::File {
        ITPFile() = default;
        ITPFile(const std::string& name) : File(name, "itp") {}
        ITPFile(const char* name) : ITPFile(std::string(name)) {}

        unsigned int size() const;

        /**
         * @brief Split a single ITP file into multiple ITP files, splitting at the specified indices. 
         * 
         * The new ITPFiles will be renumbered from 1. 
         * 
         * This is primarily intended to be used with the output from gmx::genrestr. 
         */
        std::vector<ITPFile> split_restraints(const std::vector<ITPFile>& topologies) const;
    };
}