#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace ausaxs::md {
    class make_ndx : private gmx {
        public: 
            make_ndx() {
                cmd.append("make_ndx");
            }

            make_ndx(const GROFile& gro) : make_ndx() {
                input(gro);
            }

            make_ndx& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-f", gro));
                return *this;
            }

            make_ndx& output(const NDXFile& ndx) {
                this->ndx = ndx;
                options.push_back(std::make_shared<shell::Argument>("-o", ndx));
                return *this;
            }

            std::tuple<NDXFile> run() {
                gmx::execute();
                return std::make_tuple(ndx);
            }

        private: 
            NDXFile ndx;

            void validate() const override {}

            /**
             * @brief Get the command object.
             * 
             * This is overridden since make_ndx would otherwise require user input.
             */
            shell::Command& command() override {
                return gmx::command().prepend("echo 'q' |");
            }
    };
}