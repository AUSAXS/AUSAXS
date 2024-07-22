#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <io/Folder.h>
#include <io/File.h>

namespace md {
    class genscatt : private gmx {
        public: 
            genscatt() {
                cmd.append("genscatt");
            }

            genscatt(const TPRFile& tpr, const NDXFile& ndx) : genscatt() {
                input(tpr);
                index(ndx);
            }

            genscatt& input(const TPRFile& tpr) {
                options.push_back(std::make_shared<shell::Argument>("-s", tpr));
                return *this;
            }

            genscatt& output(const std::string& prefix) {
                io::File tmp(prefix);
                this->folder = tmp.directory();
                this->prefix = tmp.stem() + tmp.extension();
                options.push_back(std::make_shared<shell::Argument>("-o", prefix));
                return *this;
            }

            genscatt& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            genscatt& group(const std::string& group) {
                cmd.prepend("echo " + group + " " + group + " |");
                return *this;
            }

            genscatt& distance(double distance) {
                options.push_back(std::make_shared<shell::Argument>("-d", distance));
                return *this;
            }

            genscatt& structure(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-s", gro));
                return *this;
            }

            std::tuple<std::vector<ITPFile>> run() {
                gmx::execute();

                std::vector<ITPFile> itps;
                for (auto const& file : folder.files()) {
                    if (file.extension() != ".itp") {continue;}
                    if (file.stem().substr(0, 6) != "scatt_") {continue;}
                    itps.emplace_back(file.path());
                }

                return std::make_tuple(itps);
            }

        private: 
            io::Folder folder;
            std::string prefix;
            void validate() const override {}
    };
}