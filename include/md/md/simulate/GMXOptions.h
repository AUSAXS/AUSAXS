#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/files/MDPCreator.h>

#include <memory>

namespace ausaxs::md {
    struct GMXOptions {
        option::Forcefield forcefield;
        option::WaterModel watermodel;
        option::BoxType boxtype;
        option::Cation cation;
        option::Anion anion;

        std::string name;
        io::Folder output;
        SHFile jobscript;
        location setupsim;
        location mainsim;
        std::shared_ptr<MDPCreator> bufmdp; // mdp file for the production simulation
        std::shared_ptr<MDPCreator> molmdp; // mdp file for the production simulation
    };

    struct BufferOptions : GMXOptions {
        BufferOptions(const GMXOptions& co, const GROFile& refgro) : GMXOptions(co), refgro(std::move(refgro)) {}
        GROFile refgro;
    };

    struct SimulateBufferOutput {
        std::unique_ptr<shell::Jobscript<MDRunResult>> job;
        TOPFile top;    // topology file
        GROFile gro;    // production structure
    };

    struct MoleculeOptions : GMXOptions {
            MoleculeOptions(const GMXOptions& co, const PDBFile& pdbfile) : GMXOptions(co), pdbfile(pdbfile) {}
            PDBFile pdbfile;
    };

    struct SimulateMoleculeOutput {
        std::unique_ptr<shell::Jobscript<MDRunResult>> job;        
        TOPFile top;    // topology file
        GROFile gro;    // solv_ion structure
    };

    struct SAXSOptions : GMXOptions {
        SAXSOptions(const GMXOptions& co, SimulateMoleculeOutput&& molecule, SimulateBufferOutput&& buffer, const PDBFile& pdb) 
            : GMXOptions(co), molecule(std::move(molecule)), buffer(std::move(buffer)), pdb(pdb) {}
        SimulateMoleculeOutput molecule;
        SimulateBufferOutput buffer;
        PDBFile pdb;
    };

    struct SAXSOutput {
        std::unique_ptr<shell::Jobscript<SAXSRunResult>> job;
    };
}