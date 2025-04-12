#pragma once

#include <md/programs/gmx.h>
#include <md/shell/Command.h>
#include <md/utility/files/all.h>
#include <md/utility/Exceptions.h>
#include <md/programs/options/forcefields/IForcefield.h>
#include <md/programs/options/water_models/IWaterModel.h>

#include <vector>

namespace ausaxs::md {
    class solvate : private gmx {
        public: 
            solvate() {
                cmd.append("solvate");
            }

            solvate(const GROFile& gro) : solvate() {
                input(gro);
            }

            solvate& input(const GROFile& gro) {
                options.push_back(std::make_shared<shell::Argument>("-cp", gro));
                return *this;
            }

            solvate& output(const GROFile& gro) {
                this->gro = gro;
                options.push_back(std::make_shared<shell::Argument>("-o", gro));
                return *this;
            }

            solvate& radius(double radius) {
                options.push_back(std::make_shared<shell::Argument>("-radius", radius));
                return *this;
            }

            solvate& topology(const TOPFile& top) {
                options.push_back(std::make_shared<shell::Argument>("-p", top));
                return *this;
            }

            solvate& solvent(observer_ptr<option::IForcefield> ff, observer_ptr<option::IWaterModel> wm) {
                wm->ensure_exists(ff);
                options.push_back(std::make_shared<shell::Argument>("-cs", wm->filename() + ".gro"));
                // options.push_back(std::make_shared<shell::Argument>("-cs", "spc216.gro"));
                return *this;
            }

            std::tuple<GROFile> run() {
                gmx::execute();
                return std::make_tuple(gro);
            }

        private: 
            GROFile gro;

            void validate() const override {}
    };
}