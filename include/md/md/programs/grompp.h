#pragma once

#include <programs/gmx.h>
#include <fstream>
#include <filesystem>
#include <utility/files/all.h>

namespace gmx {
    class grompp : private gmx {
        public: 
            grompp() {
                cmd.append("grompp");
            }

            grompp(const MDPFile& mdp, const TOPFile& top, const GROFile& gro) : grompp() {
                input(mdp);
                topology(top);
                position(gro);
            }

            grompp& input(const MDPFile& mdp) {
                options.push_back(std::make_shared<shell::Argument>("-f", mdp));
                return *this;
            }

            grompp& output(const TPRFile& tpr) {
                this->tpr = tpr;
                options.push_back(std::make_shared<shell::Argument>("-o", tpr));
                return *this;
            }

            grompp& topology(const TOPFile& top) {
                options.push_back(std::make_shared<shell::Argument>("-p", top));
                return *this;
            }

            grompp& position(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-c", gro));
                return *this;
            }

            grompp& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            grompp& restraints(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-r", gro));
                return *this;
            }

            /**
             * @brief Set the maximum number of warnings allowed.
             */
            grompp& warnings(unsigned int n) {
                options.push_back(std::make_shared<shell::Argument>("-maxwarn", n));
                return *this;
            }

            std::tuple<TPRFile> run() {
                gmx::execute();
                return std::make_tuple(std::move(tpr));
            }

        private: 
            TPRFile tpr;

            void validate() const override {}
    };
}