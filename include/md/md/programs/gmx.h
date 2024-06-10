#pragma once

#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>
#include <md/utility/files/SHFile.h>
#include <md/gmx/Settings.h>

#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <ctime>

namespace gmx {
    class gmx {
        public: 
            // shell::Command cmd = shell::Command("/data/shared/opt/gromacs/2021.5/bin/gmx");
            // shell::Command cmd = shell::Command("/data/shared/opt/gromacs/release-2021.swaxs/bin/gmx");
            shell::Command cmd = shell::Command(setting::gmx_path);

            std::vector<std::shared_ptr<shell::Option>> options;
            virtual shell::Command& command() {
                validate();
                for (auto& opt : options) {
                    cmd.append(opt->get());
                }
                cmd.append("-nocopyright -quiet");
                return cmd;
            }

            std::string execute() {
                auto cmd = command();
                if (!outputlog.empty()) {
                    cmd.prepend("set -o pipefail; ");
                    cmd.append("2>&1 | tee -a " + outputlog);                
                }
                write_cmdlog(cmd.get());

                auto result = cmd.execute();
                if (result.exit_code != 0) {
                    throw except::io_error("gmx::gmx: Error executing command: \"" + cmd.get() + "\".");
                }
                return result.out;
            }

            bool test_executable() {
                auto tmp = cmd.append("-version");
                write_cmdlog(tmp.get());
                auto res = tmp.execute();
                return res.exit_code == 0;
            };

            static void set_outputlog(const std::string& path) {
                if (std::filesystem::path(path).extension() != ".log") {
                    throw except::invalid_format("gmx::gmx: Output log file must have extension \".log\".");
                }
                std::filesystem::remove(path);
                outputlog = path;
            }

            static void set_cmdlog(const std::string& path) {
                if (std::filesystem::path(path).extension() != ".log") {
                    throw except::invalid_format("gmx::gmx: Command log file must have extension \".log\".");
                }

                cmdlog = path;
                auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                write_cmdlog(
                    "\n#################################################"
                    "\n   Program started on " + std::string(std::ctime(&time)) + 
                    "#################################################"
                );
            }

        protected:
            virtual void validate() const {}

        private: 
            inline static std::string outputlog;
            inline static std::string cmdlog;

            static void write_cmdlog(const std::string& entry) {
                if (cmdlog.empty()) {return;}
                std::ofstream log(cmdlog, std::ios_base::app);
                log << entry << std::endl;
            }

            static void write_log(const std::string& entry) {
                if (outputlog.empty()) {return;}
                std::ofstream log(outputlog, std::ios_base::app);
                log << entry << std::endl;
            }
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