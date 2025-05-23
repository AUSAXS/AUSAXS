#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/saxs.h>
#include <io/ExistingFile.h>

#include <CLI/CLI.hpp>

using namespace ausaxs;
using namespace ausaxs::md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb;
    CLI::App app{"Perform an MD buffer simulation."};
    app.add_option("input", s_pdb, "PDB structure file.")->required();
    CLI11_PARSE(app, argc, argv);

    // test executable
    if (!gmx().valid_executable()) {
        throw except::io_error("Gromacs executable not found. Please install Gromacs and add it to your PATH.");
    }

    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = s_pdb.stem(),
        .output = {"output/buffer/" + option::to_string(option::WaterModel::TIP4P2005) + "/"},
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
        .setupsim = location::lucy,
        .mainsim = location::smaug,
        // .bufmdp = std::make_shared<PRMDPCreatorSol>(),
        // .molmdp = std::make_shared<PRMDPCreatorMol>(),
        // .bufmdp = std::make_shared<FrameAnalysisMDPCreatorSol>(),
        // .molmdp = std::make_shared<FrameAnalysisMDPCreatorMol>(),
        .bufmdp = std::make_shared<TimeAnalysisMDPCreatorSol>(),
        .molmdp = std::make_shared<TimeAnalysisMDPCreatorMol>(),
    };

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(sele.output.str() + "output.log", sele.output.str() + "cmd.log");
    PDBFile pdb(s_pdb);

    GROFile conf(sele.output.str() + "setup/conf.gro");
        if (!conf.exists()) {
        // prepare the pdb file for gromacs
        auto[conf, _, posre] = pdb2gmx(pdb)
            .output({sele.output.str() + "setup/"})
            .ignore_hydrogens()
            .water_model(sele.watermodel)
            .forcefield(sele.forcefield)
        .run();
    }

    // create a box around the protein
    auto[uc] = editconf(conf)
        .output(sele.output.str() + "setup/uc.gro")
        .box_type(sele.boxtype)
        .extend(1)
    .run();

    SimulateBufferOutput buffer;
    BufferOptions bo(sele, conf);
    buffer = simulate_buffer(bo);
    buffer.job->submit();

    return 0;
} 