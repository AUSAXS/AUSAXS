#pragma once

#include <utility/files/File.h>
#include <utility/files/ITPFile.h>

#include <vector>

namespace gmx {
    // Topology file
    struct TOPFile : public detail::File {
        public:
            TOPFile() {}

            TOPFile(const std::string& name) : File(name, "top") {
                includes = discover_includes();
            }

            TOPFile(const char* name) : TOPFile(std::string(name)) {
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

            void include(const std::vector<ITPFile>& itps, const std::string& symbol);

            void fix_relative_includes();

            /**
             * @brief If the topology file contains only a single chain, extract it to a separate file.
             */
            void extract_single_chain();

            std::string copy(const Folder& folder) const override;

            const std::vector<ITPFile>& get_includes() {
                if (includes.empty()) {includes = discover_includes();}
                return includes;
            }

        private:
            std::vector<ITPFile> includes;
            std::vector<ITPFile> discover_includes() const;
            static void fix_relative_includes(const std::string& path);
    };
}