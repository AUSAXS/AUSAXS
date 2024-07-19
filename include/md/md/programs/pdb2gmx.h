#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

#include <tuple>

namespace md {
    class pdb2gmx : private gmx {
        public: 
            pdb2gmx() {
                cmd.append("pdb2gmx");
            }

            pdb2gmx(const PDBFile& pdb) : pdb2gmx() {
                input(pdb);
            }

            pdb2gmx& input(const PDBFile& pdb) {
                options.push_back(std::make_shared<shell::Argument>("-f", pdb));
                return *this;
            }

            pdb2gmx& output(const io::Folder& path) {
                folder = path;
                options.push_back(std::make_shared<shell::Argument>("-o", path + "conf.gro"));
                options.push_back(std::make_shared<shell::Argument>("-p", path + "topol.top"));
                options.push_back(std::make_shared<shell::Argument>("-i", path + "posre.itp"));
                return *this;
            }

            pdb2gmx& ignore_hydrogens() {
                options.push_back(std::make_shared<shell::Flag>("-ignh"));
                return *this;
            }

            pdb2gmx& virtual_sites() {
                options.push_back(std::make_shared<shell::Argument>("-vsite", "hydrogens"));
                return *this;
            }

            pdb2gmx& water_model(option::WaterModel) {
                // options.push_back(std::make_shared<shell::Argument>("-water", option::to_string(model)));
                cmd.prepend("echo 7 |");
                return *this;
            }

            pdb2gmx& forcefield(option::Forcefield ff) {
                options.push_back(std::make_shared<shell::Argument>("-ff", option::to_string(ff)));
                return *this;
            }

            pdb2gmx& chainsep(const std::string& sep) {
                options.push_back(std::make_shared<shell::Argument>("-chainsep", sep));
                return *this;
            }

            std::tuple<GROFile, TOPFile, ITPFile> run() {
                execute();
                TOPFile top(folder + "topol.top");
                top.extract_single_chain();
                top.fix_relative_includes();
                return std::make_tuple(GROFile(folder + "conf.gro"), top, ITPFile(folder + "posre.itp"));
            }

        private: 
            io::Folder folder;

            void validate() const override {}
    };
}