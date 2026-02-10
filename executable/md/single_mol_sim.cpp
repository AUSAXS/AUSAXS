#include <md/programs/all.h>
#include <md/simulate/SimulateBuffer.h>
#include <md/simulate/SimulateMolecule.h>
#include <md/simulate/SimulateSAXS.h>
#include <settings/MDSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/SettingsIO.h>
#include <io/ExistingFile.h>
#include <constants/Version.h>
#include <utility/Console.h>

#include <CLI/CLI.hpp>

using namespace ausaxs;
using namespace ausaxs::md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb, s_settings;
    CLI::App app{"MD simulation pipeline."};
    app.add_option("input", s_pdb, "PDB structure file.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", s_settings, "Path to the settings file.")->check(CLI::ExistingFile)->group("General options");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    CLI11_PARSE(app, argc, argv);

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(s_settings);     // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover({"."})) {
            CLI11_PARSE(app, argc, argv);
        }
    }
    if (!gmx().valid_executable()) {
        throw except::io_error("Invalid GROMACS path \"" + settings::md::gmx_path + "\"");
    }
    settings::general::output += "md/single_mol_sim/" + s_pdb.stem() + "/";

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log", settings::general::output + "cmd.log");
    PDBFile pdb(s_pdb);

    SystemSettings ss {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,
        .custom_options = {
            {"vsites", "true"},
            {"editconf extend", "0.8"},
            {"backbone_restraints", "false"},
        },
    };

    // We are essentially using GROMACS as a Monte Carlo sampler here, so we do not care about physical accuracy at all. 
    // The following settings are therefore tuned for maximum performance and fast decorrelation. 
    auto mdp = mdp::templates::production::mol()
        .add(MDPOptions::integrator = "sd")
        .add(MDPOptions::dt = "0.006")
        .add(MDPOptions::nsteps = "500000")

        // stochastic dynamics settings to get good sampling without caring about physical accuracy at all
        .add(MDPOptions::bd_fric = "10")
        .add(MDPOptions::tau_t = "0.01")
        .add(MDPOptions::ref_t = "300")
        .add(MDPOptions::tc_grps = "System")

        // constraints and vsites to allow for a large time step
        .add(MDPOptions::constraints = "all-bonds")
        .add(MDPOptions::constraint_algorithm = "lincs")
        .add(MDPOptions::lincs_order = "2")
        .add(MDPOptions::lincs_iter = "1")

        // non-bonded settings to speed up the simulation, since we do not care about accuracy
        .add(MDPOptions::coulombtype = "Reaction-Field")
        .add(MDPOptions::rcoulomb = "0.7")
        .add(MDPOptions::vdw_type = "Cutoff")
        .add(MDPOptions::rvdw = "0.7")
        .add(MDPOptions::dispcorr = "no")
        .add(MDPOptions::cutoff_scheme = "Verlet")
        .add(MDPOptions::nstlist = "100")
        .add(MDPOptions::verlet_buffer_drift = "0.01")

        // ensemble and com motion
        .add(MDPOptions::pcoupl = "no")
        .add(MDPOptions::comm_mode = "none")
        .add(MDPOptions::nstcomm = "0")

        // write frames often to get good sampling, while ensuring they are still uncorrelated
        .add(MDPOptions::nstxout_compressed = 500)
    .write(settings::general::output + "mdp/mol.mdp");
    auto molecule = simulate_molecule({
        .system = ss,
        .pdbfile = pdb,
        .mdp = mdp,
        .minimize_runner = executor::local::construct(),
        .equilibrate_runner = executor::local::construct(),
        .production_runner = executor::slurm::construct("temp/md/SmaugTemplateSampling.sh", pdb.stem() + "_mol"),
    });
    molecule.job->wait();
    return 0;
} 