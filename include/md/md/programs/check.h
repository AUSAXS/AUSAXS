#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <md/utility/Utility.h>
#include <md/utility/Exceptions.h>

namespace md {
    class check : private gmx {
        public: 
            check() {
                cmd.append("check");
            }

            check(const XTCFile& xtc) : check() {
                input(xtc);
            }

            check& input(const XTCFile& xtc) {
                options.push_back(std::make_shared<shell::Argument>("-f", xtc));
                return *this;
            }

            /**
             * @brief Determine the number of frames and the duration of a trajectory.
             * 
             * @return {number of frames, duration in ps}.
             */
            std::tuple<unsigned int, double> run() {
                auto out = gmx::execute();
                auto lines = utility::split(out, "\r\n");
                for (auto& line : lines) {
                    if (line.find("Last frame") != std::string::npos) {
                        auto tokens = utility::split(line, " \t");
                        return std::make_tuple(std::stoi(tokens[2]), std::stod(tokens[4]));
                    }
                }
                throw except::io_error("gmx::check: Could not find last frame in output.");
            }

        private: 
            void validate() const override {}
    };
}