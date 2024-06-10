#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace md {
    class select : private gmx {
        public: 
            select() {
                cmd.append("select");
            }

            select(const GROFile& gro) : select() {
                input(gro);
            }

            select& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-s", gro));
                return *this;
            }

            select& output(const NDXFile& ndx) {
                this->ndx = ndx;
                options.push_back(std::make_shared<shell::Argument>("-on", ndx));
                return *this;
            }

            select& define(const std::string& selection) {
                if (group_set) {
                    throw except::duplicate_option("select: Group already set");
                }
                group_set = true;
                cmd.prepend("echo '" + selection + "' |");
                return *this;
            }

            std::tuple<NDXFile> run() {
                gmx::execute();
                return std::make_tuple(ndx);
            }

        private: 
            bool group_set = false;
            NDXFile ndx;

            void validate() const override {
                if (!group_set) {
                    throw except::missing_option("select: Group not set");
                }
            }
    };
}