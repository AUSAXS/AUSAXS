#pragma once

#include <md/programs/all.h>
#include <md/programs/mdrun/Execution.h>
#include <md/simulate/GMXOptions.h>
#include <md/utility/files/MDPCreator.h>
#include <md/utility/Utility.h>

namespace md {
    namespace simulate {
        class Molecule {
            public:
                Molecule(MoleculeOptions& options);
                SimulateMoleculeOutput simulate();
            
            private:
                MoleculeOptions options;
                std::tuple<GROFile, TOPFile> setup();
                GROFile minimize();
                std::tuple<GROFile, NDXFile> thermalize();
                SimulateMoleculeOutput prod();
        };
    }

    SimulateMoleculeOutput simulate_molecule(MoleculeOptions& options);
}