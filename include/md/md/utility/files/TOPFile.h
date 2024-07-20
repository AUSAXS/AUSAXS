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

        void include(const std::vector<ITPFile>& itps, const std::string& symbol);

        void fix_relative_includes();

        /**
            * @brief If the topology file contains only a single chain, extract it to a separate file.
            */
        void extract_single_chain();

        std::string copy(const io::Folder& folder) const;

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