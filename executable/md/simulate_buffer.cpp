#include <md/programs/all.h>
#include <md/simulate/SimulateBuffer.h>
#include <md/simulate/SimulateMolecule.h>
#include <md/simulate/SimulateSAXS.h>
#include <io/ExistingFile.h>
#include <settings/GeneralSettings.h>

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
    settings::general::output += "md/" + s_pdb.stem() + "/";

    SystemSettings ss {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,
    };

    auto wm = option::IWaterModel::construct(ss.watermodel);
    auto ff = option::IForcefield::construct(ss.forcefield);

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log", settings::general::output + "cmd.log");
    PDBFile pdb(s_pdb);

    GROFile conf(settings::general::output + "setup/conf.gro");
        if (!conf.exists()) {
        // prepare the pdb file for gromacs
        auto[conf, _, posre] = pdb2gmx(pdb)
            .output({settings::general::output + "setup/"})
            .ignore_hydrogens()
            .water_model(wm.get())
            .forcefield(ff.get())
        .run();
    }

    // create a box around the protein
    auto[uc] = editconf(conf)
        .output(settings::general::output + "setup/uc.gro")
        .box_type(ss.boxtype)
        .extend(1)
    .run();

    auto buffer = simulate_buffer({
        .system = ss,
        .jobname = s_pdb.stem() + "_buf",
        .mdp = PRMDPCreatorSol().write(settings::general::output + "mdp/buf.mdp"),
        .setup_runner = RunLocation::local,
        .main_runner = RunLocation::local,
        .jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute_path(),
    });
    buffer.job->submit();

    return 0;
} 