#pragma once

#include <md/utility/files/ITPFile.h>
#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_top_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".top") {throw std::runtime_error("TOPFile::validate: File \"" + f->path() + "\" is not a topology file (.top).");}
            }
        };
    }

    // Topology file
    struct TOPFile : public io::detail::IValidatedFile<detail::validate_top_file> {
        using IValidatedFile::IValidatedFile;
        ~TOPFile() override = default;

        template<::detail::string_like T>
        TOPFile(const T& path) : TOPFile(std::string_view(path)) {}
        TOPFile(std::string_view name) : IValidatedFile(name) {
            includes = discover_includes();
        }

        /**
         * @brief Include an ITP file in the topology file.
         * 
         * Will be included at the end of the topology file.
         * 
         * @param itp The ITP file to include.
         * @param symbol The symbol guard. Will be used as #ifdef symbol. 
         *               If empty, no symbol guard will be used.
         */
        void include(const ITPFile& itp, const std::string& symbol);

        /**
         * @brief Include an ITP file in the topology file.
         * 
         * Will be included at the end of the specified section.
         * 
         * @param itp The ITP file to include.
         * @param symbol The symbol guard. Will be used as #ifdef symbol.
         *               If empty, no symbol guard will be used.
         * @param section The section to include the ITP file in. The first time the string is found in the file,
         *                the ITP file will be included at the end of the following block. 
         */
        void include(const ITPFile& itp, const std::string& symbol, const std::string& section);

        /**
         * @brief Include a new type of include file in this topology file. 
         *        This is useful for including e.g. backbone restraints or scattering topologies.
         *        Each chain topology include file must have a corresponding include file of the new type. 
         * 
         * @example Assuming the topology file contains the following includes:
         *          #include "topol_Protein_chain_A.itp"
         *          #include "topol_Protein_chain_B.itp"
         *
         *          Then we can add backbone restraints for each chain by calling:
         *          include_new_type({"backbone_restraints_A.itp", "backbone_restraints_B.itp"}, "POSRESBACKBONE");
         *
         *          Resulting in the following includes:
         *          #include "topol_Protein_chain_A.itp"
         *          #ifdef POSRESBACKBONE
         *              #include "backbone_restraints_A.itp"
         *          #endif
         *          ...
         *
         * @param itps The scattering ITP files to include.
         * @param symbol The symbol guard. Will be used as #ifdef symbol.
         *               If empty, no symbol guard will be used.
         */
        void include_new_type(const std::vector<ITPFile>& itps, const std::string& symbol = "");

        /**
         * @brief Fix relative includes in the topology file.
         *        Includes will be resolved relative to the topology file.
         *
         * This is useful when the topology file is moved to a different location.
         */
        void fix_relative_includes();

        /**
         * @brief If the topology file contains only a single chain, extract it to a separate file.
         *        Though this is not strictly necessary, it makes the topology file much easier to read.
         */
        void extract_single_chain();

        std::string copy(const io::Folder& folder) const;

        const std::vector<ITPFile>& get_includes() {
            if (includes.empty()) {includes = discover_includes();}
            return includes;
        }

    private:
        std::vector<ITPFile> includes;

        /**
         * @brief Scan the topology file directory for include files, and check that they are included in the topology file itself. 
         *        A warning will be printed if an include file is not included in the topology file.
         * 
         * @return A list of ITP files that are included in the topology file.
         */
        std::vector<ITPFile> discover_includes() const;

        /**
         * @brief Helper function to fix relative includes in the topology file.
         */
        static void fix_relative_includes(const io::File& path);
    };
}