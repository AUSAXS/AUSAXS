#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace ausaxs::md {
    class energy : private gmx {
        public: 
            energy() {
                cmd.append("energy");
            }

            energy(const EDRFile& edr) : energy() {
                input(edr);
            }

            energy& input(const EDRFile& edr) {
                options.push_back(std::make_shared<shell::Argument>("-f", edr));
                return *this;
            }

            energy& output(const XVGFile& xvg) {
                this->xvg = xvg;
                options.push_back(std::make_shared<shell::Argument>("-o", xvg));
                return *this;
            }

            std::tuple<XVGFile> run() {
                gmx::execute();
                return std::make_tuple(xvg);
            }

        private: 
            XVGFile xvg;
            void validate() const override {}
    };
}