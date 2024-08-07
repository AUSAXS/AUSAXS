#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <io/Folder.h>
#include <io/File.h>
#include <md/utility/Protein.h>

namespace ausaxs::md {
    class genenv : private gmx {
        public: 
            genenv() {
                cmd.append("genenv");
            }

            genenv(const XTCFile& xtc, const GROFile& gro) : genenv() {
                input(xtc);
                structure(gro);
            }

            genenv(const XTCFile& xtc, const NDXFile& ndx) : genenv() {
                input(xtc);
                index(ndx);
            }

            genenv& input(const XTCFile& xtc) {
                options.push_back(std::make_shared<shell::Argument>("-f", xtc));
                return *this;
            }

            genenv& output(const io::Folder& folder) {
                this->folder = folder;
                return *this;
            }

            genenv& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            genenv& distance(double distance) {
                options.push_back(std::make_shared<shell::Argument>("-d", distance));
                return *this;
            }

            genenv& structure(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-s", gro));
                return *this;
            }

            genenv& group(const std::string& group) {
                cmd.prepend("echo " + group + " " + group + " |");
                return *this;
            }

            /**
             * @brief Generate a spherical envelope around the structure.
             *        Note that when using this option, the input trajectory and structure are ignored. 
             */
            genenv& sphere(double radius) {
                options.push_back(std::make_shared<shell::Flag>("-sphere"));
                options.push_back(std::make_shared<shell::Argument>("-d_sphere", radius));
                return *this;
            }

            /**
             * @brief Generate a spherical envelope around the structure.
             *        Note that when using this option, the input trajectory and structure are ignored. 
             *
             * @param gro: The structure file. The radius will be added to the half of the maximum extent of the structure.
             * @param radius: The additional extent of the radius of the sphere.
             */
            genenv& sphere(const PDBFile& gro, double radius) {
                auto dmax = Protein(gro).maximum_distance()/20 + 0.2; // half of the maximum extent converted to nm plus generous 2Å vdw buffer 
                options.push_back(std::make_shared<shell::Flag>("-sphere"));
                options.push_back(std::make_shared<shell::Argument>("-d_sphere", dmax + radius));
                return *this;
            }

            /**
             * @return unsigned int: The good pbc index
             *         GROFile: The envelope gro file
             *         PYFile: The envelope python file
             *         DATFile: The envelope dat file
             */
            std::tuple<unsigned int, GROFile, PYFile, DATFile> run() {
                auto out = gmx::execute();

                GROFile gro;
                PYFile py;
                DATFile dat;
                io::Folder root(".");
                for (auto& entry : root.files()) {
                    auto path = entry.path();
                    if (path.find("envelope") != std::string::npos) {
                        if (path.find(".gro") != std::string::npos) {
                            gro = GROFile(path);
                        } else if (path.find(".py") != std::string::npos) {
                            py = PYFile(path);
                        } else if (path.find(".dat") != std::string::npos) {
                            dat = DATFile(path);
                        }
                    }
                }

                // move to output folder
                if (!folder.path().empty()) {
                    if (gro.exists()) {gro.move(folder);}
                    py.move(folder);
                    dat.move(folder);
                }

                // find good pbc in output
                std::string pbc;
                unsigned int pbc_start = out.find("Global atom number = ") + 21;
                for (unsigned int i = pbc_start; i < pbc_start+10; i++) {
                    if (out[i] == ' ') {continue;}
                    if (std::isdigit(out[i])) {pbc += out[i];}
                    else {break;}
                }

                if (pbc.empty()) {return std::make_tuple(0, gro, py, dat);}
                return std::make_tuple(std::stoi(pbc), gro, py, dat);
            }

        private: 
            io::Folder folder;
            void validate() const override {}
    };
}