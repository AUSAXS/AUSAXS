#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace ausaxs::md {
    class editconf : private gmx {
        public:
            editconf() {
                cmd.append("editconf");
            }

            editconf(const GROFile& gro) : editconf() {
                input(gro);
            }

            editconf& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-f", gro));
                return *this;
            }

            editconf& output(const GROFile& gro) {
                this->gro = gro;
                options.push_back(std::make_shared<shell::Argument>("-o", gro));
                return *this;
            }

            editconf& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            editconf& select(const std::string& group) {
                options.push_back(std::make_shared<shell::Argument>("-ndef", group));
                return *this;
            }

            editconf& box_type(option::BoxType type) {
                options.push_back(std::make_shared<shell::Argument>("-bt", option::to_string(type)));
                return *this;
            }

            editconf& extend(double dist) {
                options.push_back(std::make_shared<shell::Argument>("-d", dist));
                return *this;
            }

            std::tuple<GROFile> run() {
                gmx::execute();
                return std::make_tuple(std::move(gro));
            }

        private: 
            GROFile gro;

            void validate() const override {}
    };
}