#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/timeanalysis.h>
#include <io/ExistingFile.h>
#include <settings/GeneralSettings.h>

#include <CLI/CLI.hpp>

using namespace ausaxs;
using namespace ausaxs::md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb;
    CLI::App app{"Perform an analysis of the SAXS-dependency on the simulation time."};
    app.add_option("input", s_pdb, "PDB structure file.")->required();
    CLI11_PARSE(app, argc, argv);

    // test executable
    if (!gmx().valid_executable()) {
        throw except::io_error("Gromacs executable not found. Please install Gromacs and add it to your PATH.");
    }
    settings::general::output += "md/timeanalysis/" + s_pdb.stem() + "/";

    SystemSettings ss {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,
    };

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log", settings::general::output + "cmd.log");
    PDBFile pdb(s_pdb);

    auto molecule = simulate_molecule({
        .system = ss,
        .jobname = s_pdb.stem() + "_mol",
        .pdbfile = pdb,
        .mdp = TimeAnalysisMDPCreatorMol().write(settings::general::output + "mdp/ta_mol.mdp"),
        .setup_runner = RunLocation::local,
        .main_runner = RunLocation::local,
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
    });

    auto buffer = simulate_buffer({
        .system = ss,
        .jobname = s_pdb.stem() + "_buf",
        .refgro = molecule.gro,
        .mdp = TimeAnalysisMDPCreatorSol().write(settings::general::output + "mdp/ta_buf.mdp"),
        .setup_runner = RunLocation::local,
        .main_runner = RunLocation::local,
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
    });

    auto saxs = timeanalysis({
        .pdbfile = pdb,
        .molecule = std::move(molecule),
        .buffer = std::move(buffer),
        .runner = RunLocation::local,
        .jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute_path(),
    }, 500);

    for (int i = 0; i < static_cast<int>(saxs.size()); i++) {saxs[i].job->submit();}
    for (int i = static_cast<int>(saxs.size())-1; i >= 0; i--) { // reverse wait order since the last job is the shortest run
        saxs[i].job->wait();
        auto res = saxs[i].job->result();
        res.xvg.rename(std::to_string(i+1) + ".xvg").copy({settings::general::output + "output"});
    }
    return 0;
} 