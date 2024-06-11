#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace md {
    class genion : private gmx {
        public: 
            genion() {
                cmd.append("genion");
            }

            genion(const TPRFile& tpr) : genion() {
                input(tpr);
            }

            genion& input(const TPRFile& tpr) {
                options.push_back(std::make_shared<shell::Argument>("-s", tpr));
                return *this;
            }

            genion& output(const GROFile& gro) {
                this->gro = gro;
                options.push_back(std::make_shared<shell::Argument>("-o", gro));
                return *this;
            }

            genion& topology(const TOPFile& top) {
                options.push_back(std::make_shared<shell::Argument>("-p", top));
                return *this;
            }

            genion& neutralize() {
                options.push_back(std::make_shared<shell::Flag>("-neutral"));
                return *this;
            }

            genion& anion(option::Anion anion) {
                options.push_back(std::make_shared<shell::Argument>("-nname", option::to_string(anion)));
                return *this;
            }

            genion& cation(option::Cation cation) {
                options.push_back(std::make_shared<shell::Argument>("-pname", option::to_string(cation)));
                return *this;
            }

            /**
             * @brief Set the index group to add the ions to.
             */
            genion& ion_group(std::string group) {
                cmd.prepend("echo \"" + group + "\" |");
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