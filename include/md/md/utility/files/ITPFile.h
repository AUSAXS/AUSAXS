#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_itp_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".itp") {throw std::runtime_error("ITPFile::validate: File \"" + f->path() + "\" is not a include topology file (.itp).");}
            }
        };
    }

    // Include topology file
    struct ITPFile : public io::detail::IValidatedFile<detail::validate_itp_file> {
        using IValidatedFile::IValidatedFile;

        unsigned int size() const;

        /**
         * @brief Split a single ITP file into multiple ITP files, splitting at the specified indices. 
         * 
         * The new ITPFiles will be renumbered from 1. 
         * This is primarily intended to be used with the output from gmx::genrestr. 
         */
        std::vector<ITPFile> split_restraints(const std::vector<ITPFile>& topologies) const;

        /**
         * @brief Standardize the name of this and all nested ITP files.
         * 
         * The name will be changed to e.g. "topol_prefix" or "backbone_prefix".
         */
        void standardize_name(std::string_view postfix);
    };
}