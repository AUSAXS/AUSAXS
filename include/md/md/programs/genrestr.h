#pragma once

#include <programs/gmx.h>
#include <utility/files/all.h>

namespace gmx {
    class genrestr : private gmx {
        public:
            genrestr() {
                cmd.append("genrestr");
            }

            genrestr(const GROFile& gro) : genrestr() {
                input(gro);
            }

            genrestr& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-f", gro));
                return *this;
            }

            genrestr& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            genrestr& output(const ITPFile& itp) {
                this->itp = itp;
                options.push_back(std::make_shared<shell::Argument>("-o", itp));
                return *this;
            }

            genrestr& force(int x, int y, int z) {
                options.push_back(std::make_shared<shell::Argument>("-fc", std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z)));
                return *this;
            }

            /**
             * @brief Get the command object.
             * 
             * This is overridden since genrestr would otherwise require user input.
             * Currently only backbone restraints are supported.
             */
            shell::Command& command() override {
                return gmx::command().prepend("echo 'Backbone' |");
            }

            std::tuple<ITPFile> run() {
                gmx::execute();
                return std::make_tuple(std::move(itp));
            }

        private: 
            ITPFile itp;

            void validate() const override {}
    };
}