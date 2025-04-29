#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <md/programs/options/water_models/IWaterModel.h>
#include <md/programs/options/forcefields/IForcefield.h>

#include <tuple>

namespace ausaxs::md {
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
                options.push_back(std::make_shared<shell::Argument>("-o", path.str() + "/conf.gro"));
                options.push_back(std::make_shared<shell::Argument>("-p", path.str() + "/topol.top"));
                options.push_back(std::make_shared<shell::Argument>("-i", path.str() + "/posre.itp"));
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

            pdb2gmx& water_model(observer_ptr<option::IWaterModel> wm) {
                this->wm = wm;
                return *this;
            }

            pdb2gmx& forcefield(observer_ptr<option::IForcefield> ff) {
                this->ff = ff;
                options.push_back(std::make_shared<shell::Argument>("-ff", ff->filename()));
                return *this;
            }

            pdb2gmx& chainsep(const std::string& sep) {
                options.push_back(std::make_shared<shell::Argument>("-chainsep", sep));
                return *this;
            }

            std::tuple<GROFile, TOPFile, ITPFile> run() {
                // prepend water selection
                wm->ensure_exists(ff);
                {
                    auto res = shell::Command("cat " + settings::md::gmx_top_path + ff->filename() + ".ff/watermodels.dat").execute();
                    if (res.exit_code != 0) {
                        throw except::io_error("pdb2gmx: Error executing command: \"" + res.out + "\".");
                    }
                    std::string wm_s = wm->filename();
    
                    std::vector<std::string> lines;
                    int line_start = 0, line_end = 0;
                    for (auto& c : res.out) {
                        if (c == '\n') {
                            lines.push_back(res.out.substr(line_start, line_end - line_start));
                            ++line_end;
                            line_start = line_end;
                        } else {
                            ++line_end;
                        }
                    }
                    lines.push_back(res.out.substr(line_start, line_end - line_start));
                    int c = 1;
                    bool found = false;
                    for (auto& line : lines) {
                        if (line.starts_with(wm_s)) {
                            cmd.prepend("echo " + std::to_string(c) + " |");
                            found = true;
                            break;
                        }
                        ++c;
                    }
                    if (!found) {
                        throw except::io_error(
                            "pdb2gmx: Could not find water model index \"" + wm_s + "\" in file \"" 
                            + settings::md::gmx_top_path + ff->filename() + ".ff/watermodels.dat\"."
                        );
                    }
                }

                execute();
                TOPFile top(folder.str() + "/topol.top");
                top.extract_single_chain();
                top.fix_relative_includes();
                top.standardize_itp_names();
                return std::make_tuple(GROFile(folder.str() + "/conf.gro"), top, ITPFile(folder.str() + "/posre.itp"));
            }

        private:
            io::Folder folder;
            observer_ptr<option::IForcefield> ff;
            observer_ptr<option::IWaterModel> wm;

            void validate() const override {}
    };
}