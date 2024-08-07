#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/saxs.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <constants/Version.h>
#include <utility/Console.h>

#include <CLI/CLI.hpp>

using namespace md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb, s_settings;
    CLI::App app{"Extract a water profile from GROMACS."};
    app.add_option("input", s_pdb, "PDB structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--buffer", settings::md::buffer_path, "Pre-simulated buffer. Set this to the base output folder.");
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

    settings::general::output = "output/md/" + s_pdb.stem() + "/water_profile/";

    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = s_pdb.stem(),
        .output = {"output/md/" + s_pdb.stem() + "/"},
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
        .setupsim = location::local,
        .mainsim = location::local,
        .bufmdp = std::make_shared<PRMDPCreatorSol>(),
        .molmdp = std::make_shared<PRMDPCreatorMol>(),
    };

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log", settings::general::output + "cmd.log");
    PDBFile pdb(s_pdb);

    //##################################//
    //###           SETUP            ###//
    //##################################//
    console::print_info("\nPreparing simulation for " + s_pdb.filename());
    console::indent();

    io::Folder mdp_folder = settings::general::output + "mdp/";
    io::Folder setup_path = settings::general::output + "protein/setup/";
    io::Folder em_path    = settings::general::output + "protein/em/";
    io::Folder eq_path    = settings::general::output + "protein/eq/";
    io::Folder prod_path  = settings::general::output + "protein/prod/";
    mdp_folder.create(); setup_path.create(); em_path.create(); eq_path.create(); prod_path.create();

    GROFile solv_ion(setup_path + "solv_ion.gro");
    TOPFile top(setup_path + "topol.top");
    if (!solv_ion.exists() || !top.exists()) {
        GROFile conf(setup_path + "conf.gro");
        if (!conf.exists()) {
            console::print_text("Converting structure to GROMACS format...");

            // prepare the pdb file for gromacs
            auto[conf, top, posre] = pdb2gmx(pdb)
                .output(setup_path)
                .ignore_hydrogens()
                // .virtual_sites()
                .water_model(option::WaterModel::TIP4P2005)
                .forcefield(option::Forcefield::AMBER99SB_ILDN)
            .run();
        } else {
            console::print_text("Reusing previously generated GROMACS structure.");
        }

        // create a box around the protein
        console::print_text("Generating unit cell...");
        auto[uc] = editconf(conf)
            .output(setup_path + "uc.gro")
            .box_type(option::BoxType::DODECAHEDRON)
            .extend(1.5)
        .run();

        // add water to the box
        console::print_text("Solvating unit cell...");
        GROFile solv;
        {
            // dummy solvent purely to get the number of solvent molecules
            auto[tmp] = solvate(uc)
                .output(setup_path + "solv.gro")
                .solvent(option::Forcefield::AMBER99SB_ILDN, option::WaterModel::TIP4P2005)
                .radius(0.2)
                .topology(top)
            .run();

            // generate the actual solvent in random positions
            std::tie(solv) = insert_molecules(uc)
                .solvent(option::Forcefield::AMBER99SB_ILDN, option::WaterModel::TIP4P2005)
                .nmol(tmp.size_solvent())
            .run();
        }

        // generate an empty tpr file 
        MDPFile mdp = MDPFile(setup_path + "empty.mdp"); mdp.create();
        auto[ions] = grompp(mdp, top, solv)
            .output(setup_path + "ions.tpr")
            .warnings(1)
        .run();
        mdp.remove();

        // add ions to the box
        console::print_text("Adding ions to unit cell...");
        genion(ions)
            .output(solv_ion)
            .topology(top)
            .neutralize()
            .cation(option::Cation::NA)
            .anion(option::Anion::CL)
            .ion_group("SOL")
        .run();
    } else {
        console::print_text("Reusing previously generated system setup.");
    }

    return 0;
}