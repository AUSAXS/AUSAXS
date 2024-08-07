#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <md/utility/Exceptions.h>
#include <utility/StringUtils.h>

namespace md {
    class insert_molecules : private gmx {
        public: 
            insert_molecules() {
                cmd.append("insert-molecules");
            }

            insert_molecules(const GROFile& gro) : insert_molecules() {
                input(gro);
            }

            insert_molecules& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-f", gro));
                return *this;
            }

            insert_molecules& solvent(const option::Forcefield ff, const option::WaterModel& solv) {
                options.push_back(std::make_shared<shell::Argument>("-ci", option::to_string(ff) + ".ff/" + option::to_string(solv) + ".gro"));
                return *this;
            }

            insert_molecules& output(const io::File& gro) {
                options.push_back(std::make_shared<shell::Argument>("-o", gro));
                out = gro;
                return *this;
            }

            insert_molecules& nmol(unsigned int n) {
                options.push_back(std::make_shared<shell::Argument>("-nmol", n));
                return *this;
            }

            std::tuple<GROFile> run() {
                gmx::execute();
                return std::make_tuple(out);
            }

        private: 
            io::File out;
            void validate() const override {}
    };
}