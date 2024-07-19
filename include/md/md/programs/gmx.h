#pragma once

#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>
#include <md/utility/files/SHFile.h>
#include <settings/MDSettings.h>

#include <string>
#include <vector>

namespace md {
    class gmx {
        public: 
            shell::Command cmd = shell::Command(settings::md::gmx_path);

            std::vector<std::shared_ptr<shell::Option>> options;
            virtual shell::Command& command();

            std::string execute();

            bool valid_executable();

            static void set_outputlog(const std::string& path);
            static void set_cmdlog(const std::string& path);

        protected:
            virtual void validate() const;

        private: 
            inline static std::string outputlog;
            inline static std::string cmdlog;

            static void write_cmdlog(const std::string& entry);
            static void write_log(const std::string& entry);
    };

    namespace option {
        enum class Forcefield {
            AMBER99SB,
            AMBER99SB_ILDN
        };

        enum class WaterModel {
            TIP3P,
            TIP4P,
            TIP4P2005
        };

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

        std::string to_string(Forcefield opt);
        std::string to_string(WaterModel opt);
        std::string to_string(BoxType opt);
        std::string to_string(Cation opt);
        std::string to_string(Anion opt);
    }
}