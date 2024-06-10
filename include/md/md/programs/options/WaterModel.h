#pragma once

#include <md/programs/all.h>
#include <md/utility/files/MDPCreator.h>

namespace gmx::option {
    struct WaterModel {
        virtual ~WaterModel() = default;

        virtual std::string name() const = 0;
        virtual void generate_mdp(const GROFile& gro) const = 0;
    };

    struct TIP4P2005 : WaterModel {
        std::string name() const override { return "tip4p2005"; }
        void generate_mdp(const GROFile& gro) const override {
            // generate protein index file
            auto[tmp1] = select(gro)
                .output("tmp.ndx")
                .define("'OnlyWater' name 'OW' or name 'HW1' or name 'HW2' or name 'HW3' or group 'Ion'")
            .run();

            auto[molind] = make_ndx(gro)
                .output("protein/index.ndx")
            .run();

            molind.append_file(tmp1);
            tmp1.remove();

            // generate buffer index file
            auto[tmp2] = select(gro)
                .output("tmp.ndx")
                .define("'OnlyWater' name 'OW' or name 'HW1' or name 'HW2' or name 'HW3'")
            .run();

            auto[bufind] = make_ndx(gro)
                .output("buffer/index.ndx")
            .run();

            bufind.append_file(tmp2);
            tmp2.remove();

            // generate mdp file
            auto molmdp = SAXSMDPCreatorMol();
            auto solmdp = SAXSMDPCreatorSol();
        }
    };
}