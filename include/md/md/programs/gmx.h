#pragma once

#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>
#include <md/utility/files/SHFile.h>
#include <settings/MDSettings.h>

#include <string>
#include <vector>

namespace ausaxs::md {
    class gmx {
        public: 
            shell::Command cmd = shell::Command(settings::md::gmx_path);

            std::vector<std::shared_ptr<shell::Option>> options;
            virtual shell::Command& command();

            std::string execute();

            bool valid_executable();

            static void set_logfile(const io::File& log, const io::File& cmdlog);

        protected:
            virtual void validate() const;

        private: 
            inline static bool log = false;
            inline static io::File outputlog;
            inline static io::File cmdlog;

            static void write_cmdlog(std::string_view entry);
            static void write_log(std::string_view entry);
    };

    namespace option {
        enum class BoxType {
            CUBIC,
            TRICLINIC,
            DODECAHEDRON,
            OCTAHEDRON
        };

        enum class Cation {
            NA
        };
        
        enum class Anion {
            CL
        };

        std::string to_string(BoxType opt);
        std::string to_string(Cation opt);
        std::string to_string(Anion opt);
    }
}