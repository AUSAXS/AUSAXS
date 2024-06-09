#pragma once

#include <programs/gmx.h>
#include <utility/files/all.h>

#include <iostream>

namespace gmx {
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
                this->folder = detail::File(prefix).parent_path();
                this->prefix = detail::File(prefix).filename();
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
                for (auto const& entry : std::filesystem::directory_iterator(folder.path)) {
                    auto path = entry.path();
                    if (path.extension() != ".itp") {continue;}
                    if (path.stem().string().substr(0, 6) != "scatt_") {continue;}
                    itps.emplace_back(path.string());
                }

                return std::make_tuple(itps);
            }

        private: 
            Folder folder;
            std::string prefix;
            void validate() const override {}
    };
}