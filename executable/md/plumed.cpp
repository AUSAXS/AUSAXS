#include <md/simulate/GMXOptions.h>
#include <md/programs/all.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <settings/MDSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/SettingsIO.h>
#include <io/ExistingFile.h>
#include <constants/Version.h>
#include <utility/Console.h>

#include <CLI/CLI.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::md;

static std::pair<int, int> find_terminal_ca_atoms(const GROFile& gro) {
    std::ifstream f(gro.absolute_path());
    if (!f) {throw std::runtime_error("find_terminal_ca_atoms: cannot open " + gro.path());}

    std::string line;
    std::getline(f, line); // title line
    std::getline(f, line); // atom count line

    int first_ca = -1, last_ca = -1, atom_idx = 0;
    while (std::getline(f, line)) {
        // The last (box) line is shorter than a real atom line; stop there.
        if (line.size() < 20) {break;}
        ++atom_idx;

        // GRO format: cols 10-14 (5 chars) hold the atom name, left-padded.
        std::string atom_name = line.substr(10, 5);
        // Trim whitespace
        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(' '));
            auto last = s.find_last_not_of(' ');
            if (last != std::string::npos) {s.erase(last + 1);}
        };
        trim(atom_name);

        if (atom_name == "CA") {
            if (first_ca == -1) {first_ca = atom_idx;}
            last_ca = atom_idx;
        }
    }
    return {first_ca, last_ca};
}

// ---------------------------------------------------------------------------
// Generate the PLUMED input file content.
//
// The bias computes the end-to-end distance (atom1 to atom2) for each replica,
// averages it with ENSEMBLE across all N replicas, then applies a MaxEnt
// restraint that drives <ete> → target_d_nm.
//
// MAXENT parameter guide:
//   AT     - target value for the ensemble-averaged CV (nm)
//   KAPPA  - learning rate for the Lagrange multiplier λ (kJ/mol/nm^2 × ps^-1)
//            A larger value tracks the target more aggressively.
//   TEMP   - simulation temperature (K)
//   TAU    - memory time for the running average of <CV> (ps); controls the
//            timescale over which ensemble averages are accumulated.
// ---------------------------------------------------------------------------
static std::string make_plumed_input(
        int atom1, int atom2,
        double target_d_nm, double kappa, double temp,
        int n_replicas)
{
    std::ostringstream s;
    s << "# PLUMED input: MaxEnt bias on the ensemble-averaged end-to-end distance\n"
      << "# N replicas: " << n_replicas << "  |  target <ete>: " << target_d_nm << " nm\n\n";

    // Keep molecule whole across PBC up to the last end-to-end atom.
    s << "WHOLEMOLECULES ENTITY0=1-" << atom2 << "\n\n";

    // Per-replica end-to-end distance (Cα N-term → Cα C-term).
    s << "# End-to-end distance for this replica\n"
      << "ete: DISTANCE ATOMS=" << atom1 << "," << atom2 << "\n\n";

    // Replica-averaged CV.
    s << "# Ensemble average over all " << n_replicas << " replicas\n"
      << "ens_ete: ENSEMBLE ARG=ete\n\n";

    // Maximum entropy restraint.
    s << "# MaxEnt restraint: adapts lambda so that <ete> = " << target_d_nm << " nm\n"
      << "MAXENT ...\n"
      << "  LABEL=maxent\n"
      << "  TYPE=GAUSS\n"
      << "  ARG=ens_ete.ete\n"
      << "  AT=" << target_d_nm << "\n"
      << "  KAPPA=" << kappa << "\n"
      << "  TEMP=" << temp << "\n"
      << "  TAU=1000.0\n"
      << "... MAXENT\n\n";

    // Output.
    s << "# Write per-replica distance, ensemble average, and lambda to COLVAR\n"
      << "PRINT ARG=ete,ens_ete.ete,maxent.* STRIDE=100 FILE=COLVAR\n";

    return s.str();
}

// ===========================================================================
int main(int argc, char const* argv[]) {
    io::ExistingFile s_pdb, s_settings;
    int    N           = 4;
    long   nsteps      = 500000;  // prod steps; override with --nsteps for quick debug runs
    double target_d    = 3.0;   // nm
    double kappa       = 10.0;  // kJ/mol/nm² per ps (lambda learning rate)
    double temperature = 300.0; // K

    CLI::App app{"GROMACS+PLUMED multi-replica MaxEnt end-to-end distance bias."};
    app.add_option("input", s_pdb, "PDB structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("-N,--replicas",  N,         "Number of replicas (default: 4).")->default_val(4);
    app.add_option("-d,--distance",  target_d,  "Target ensemble-averaged end-to-end distance (nm).")->required();
    app.add_option("--nsteps",       nsteps,    "Production steps (default: 500000; use a small value for debugging).");
    app.add_option("--kappa",        kappa,     "MaxEnt lambda learning rate (kJ/mol/nm²/ps, default: 10).")->default_val(10.0);
    app.add_option("--temp",         temperature, "Simulation temperature in K (default: 300).")->default_val(300.0);

    auto p_settings = app.add_option("-s,--settings", s_settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    CLI11_PARSE(app, argc, argv);

    if (p_settings->count() != 0) {
        settings::read(s_settings);
        CLI11_PARSE(app, argc, argv);
    } else {
        if (settings::discover({"."})) {
            CLI11_PARSE(app, argc, argv);
        }
    }

    if (!gmx().valid_executable()) {
        throw except::io_error("Invalid GROMACS path \"" + settings::md::gmx_path + "\"");
    }

    settings::general::output += "md/plumed_maxent/" + s_pdb.stem() + "_" + std::to_string(N) + "/";
    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log",
                           settings::general::output + "cmd.log");

    PDBFile pdb(s_pdb);
    console::print_info("\nMaxEnt multi-replica simulation for " + pdb.filename() + " (" + std::to_string(N) + " replicas, target d = " + std::to_string(target_d) + " nm)");

    std::unique_ptr<executor::type> em_runner = executor::templated::construct("temp/md/lusiTemplate.sh");
    std::unique_ptr<executor::type> eq_runner = executor::templated::construct("temp/md/lusiTemplate.sh");
    std::unique_ptr<executor::type> prod_runner = executor::multi::slurm(N, "temp/md/SmaugTemplateMPI.sh", "plumed");

    //##################################//
    //###      SYSTEM PARAMETERS     ###//
    //##################################//
    SystemSettings ss {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype    = option::BoxType::DODECAHEDRON,
        .cation     = option::Cation::NA,
        .anion      = option::Anion::CL,
        .custom_options = {
            {"vsites",            "true"},
            {"editconf extend",   "0.85"},
            {"backbone_restraints","false"},
        },
    };

    auto ff = option::IForcefield::construct(ss.forcefield);
    auto wm = option::IWaterModel::construct(ss.watermodel);

    //##################################//
    //###        PRODUCTION MDP      ###//
    //##################################//
    // Stochastic dynamics (sd) for good conformational sampling without
    // caring about physical accuracy; large time step via virtual sites and
    // all-bonds constraints. See single_mol_sim.cpp for rationale.
    auto mdp_prod = mdp::templates::production::mol()
        .add(MDPOptions::integrator          = "sd")
        .add(MDPOptions::dt                  = "0.004")
        .add(MDPOptions::nsteps              = std::to_string(nsteps))
        .add(MDPOptions::define              = "")   // remove posresbackbone

        .add(MDPOptions::tau_t               = "3")
        .add(MDPOptions::ref_t               = "300")
        .add(MDPOptions::tc_grps             = "System")

        .add(MDPOptions::constraints         = "all-bonds")
        .add(MDPOptions::constraint_algorithm= "lincs")
        .add(MDPOptions::lincs_order         = "4")
        .add(MDPOptions::lincs_iter          = "2")

        .add(MDPOptions::coulombtype         = "Reaction-Field")
        .add(MDPOptions::rcoulomb            = "0.8")
        .add(MDPOptions::vdw_type            = "Cutoff")
        .add(MDPOptions::rvdw                = "0.8")
        .add(MDPOptions::dispcorr            = "no")
        .add(MDPOptions::cutoff_scheme       = "Verlet")
        .add(MDPOptions::nstlist             = "100")
        .add(MDPOptions::verlet_buffer_drift = "0.01")

        .add(MDPOptions::pcoupl              = "no")
        .add(MDPOptions::tcoupl              = "no")
        .add(MDPOptions::comm_mode           = "linear")
        .add(MDPOptions::nstcomm             = "100")

        .add(MDPOptions::nstxout_compressed  = 5000)
    ;

    //##################################//
    //###           SETUP            ###//
    //##################################//
    io::Folder setup_path = settings::general::output + "setup/";
    io::Folder em_path    = settings::general::output + "em/";
    io::Folder eq_path    = settings::general::output + "eq/";
    io::Folder mdp_folder = settings::general::output + "mdp/";
    setup_path.create(); em_path.create(); eq_path.create(); mdp_folder.create();

    GROFile solv_ion(setup_path + "solv_ion.gro");
    TOPFile top(setup_path + "topol.top");

    if (!solv_ion.exists() || !top.exists()) {
        GROFile conf(setup_path + "conf.gro");
        if (!conf.exists()) {
            console::print_text("Converting structure to GROMACS format...");
            auto cmd = pdb2gmx(pdb)
                .output(setup_path)
                .ignore_hydrogens()
                .water_model(wm.get())
                .forcefield(ff.get())
                .virtual_sites();
            auto [conf_out, top_out, posre] = cmd.run();
        } else {
            console::print_text("Reusing previously generated GROMACS structure.");
        }

        console::print_text("Generating unit cell...");
        auto [uc] = editconf(conf)
            .output(setup_path + "uc.gro")
            .box_type(ss.boxtype)
            .extend(0.85)
        .run();

        console::print_text("Solvating unit cell...");
        auto [solv] = solvate(uc)
            .output(setup_path + "solv.gro")
            .solvent(ff.get(), wm.get())
            .radius(0.2)
            .topology(top)
        .run();

        // Temporary empty MDP for ion addition
        MDPFile mdp_empty(setup_path + "empty.mdp"); mdp_empty.create();
        auto [ions_tpr] = grompp(mdp_empty, top, solv)
            .output(setup_path + "ions.tpr")
            .warnings(1)
        .run();
        mdp_empty.remove();

        console::print_text("Adding ions...");
        genion(ions_tpr)
            .output(solv_ion)
            .topology(top)
            .neutralize()
            .cation(ss.cation)
            .anion(ss.anion)
            .ion_group("SOL")
        .run();
    } else {
        console::print_text("Reusing previously generated system setup.");
    }

    //##################################//
    //###     ENERGY MINIMIZATION    ###//
    //##################################//
    GROFile emgro(em_path + "em.gro");
    if (!emgro.exists()) {
        console::print_text("Running energy minimization...");
        MDPFile mdp_em = mdp::templates::minimize::base().write(mdp_folder + "em.mdp");
        auto [emtpr] = grompp(mdp_em, top, solv_ion)
            .output(em_path + "em.tpr")
            .restraints(solv_ion)
        .run();
        mdrun(emtpr)
            .output(em_path, "/em")
        .run(std::move(em_runner))->wait();
    } else {
        console::print_text("Reusing previously generated energy minimization.");
    }

    //##################################//
    //###       EQUILIBRATION        ###//
    //##################################//
    GROFile eqgro(eq_path + "eq.gro");
    NDXFile index(eq_path + "index.ndx");
    if (!eqgro.exists()) {
        console::print_text("Running equilibration...");
        auto [_] = make_ndx(emgro).output(index).run();

        // Ensure "Water_and_Ions" group exists (may be absent if no ions were added).
        if (!index.contains("Water_and_Ions")) {
            auto [tmp] = md::select(emgro)
                .output("tmp.ndx")
                .define("\"Water_and_Ions\" group \"Water\"")
            .run();
            index.append_file(tmp);
            tmp.remove();
        }

        MDPFile mdp_eq = mdp_prod
            .add(MDPOptions::nsteps = "50000")
        .write(mdp_folder + "eq.mdp");

        auto [eqtpr] = grompp(mdp_eq, top, emgro)
            .output(eq_path + "eq.tpr")
            .index(index)
            .restraints(emgro)
            .warnings(1)
        .run();
        mdrun(eqtpr)
            .output(eq_path, "/eq")
        .run(std::move(eq_runner))->wait();
    } else {
        console::print_text("Reusing previously generated equilibration.");
    }

    //##################################//
    //###   END-TO-END ATOM INDICES  ###//
    //##################################//
    auto [ca_first, ca_last] = find_terminal_ca_atoms(eqgro);
    if (ca_first < 0 || ca_last < 0) {
        throw std::runtime_error(
            "Could not identify Cα atoms for end-to-end distance. "
            "Please specify --atom1 and --atom2 explicitly.");
    }
    console::print_text("End-to-end atoms: index " + std::to_string(ca_first)
        + " (first Cα) and " + std::to_string(ca_last) + " (last Cα)");

    //##################################//
    //###        PLUMED FILE         ###//
    //##################################//
    io::File plumed_file(settings::general::output + "plumed.dat");
    {
        std::ofstream out(plumed_file.path());
        if (!out) {throw std::runtime_error("Could not create PLUMED input file at " + plumed_file.path());}
        out << make_plumed_input(ca_first, ca_last, target_d, kappa, temperature, N);
    }
    console::print_text("PLUMED input written to " + plumed_file.path());

    //##################################//
    //###   REPLICA DIRECTORIES      ###//
    //##################################//
    // Each replica gets its own sub-directory containing a TPR file.
    // GROMACS -multidir requires all TPR files to have the same deffnm.
    io::Folder prod_base(settings::general::output + "prod/");
    auto mdp_prod_file = mdp_prod.write(mdp_folder + "prod.mdp");
    prod_base.create();

    std::vector<io::Folder> rep_dirs;

    for (int i = 0; i < N; ++i) {
        io::Folder rep_dir(prod_base + ("rep" + std::to_string(i)));
        rep_dir.create();

        TPRFile rep_tpr(rep_dir + "prod.tpr");
        if (!rep_tpr.exists()) {
            grompp(mdp_prod_file, top, eqgro)
                .output(rep_tpr)
                .index(index)
                .warnings(1)
            .run();
        }

        rep_dirs.push_back(rep_dir);
    }

    //##################################//
    //###   MULTI-REPLICA MDRUN      ###//
    //##################################//
    // gmx mdrun with -multidir sets up N independent replicas communicating
    // through PLUMED. Each replica reads prod.tpr from its own directory and
    // writes output files (prod.xtc, prod.gro, …) there.
    //
    // The PLUMED ENSEMBLE action exchanges CV values between replicas every
    // step, so they must all be on the same compute node (MPI ranks) or
    // connected via the PLUMED multi-sim infrastructure.
    console::print_info("Launching " + std::to_string(N)
        + "-replica PLUMED MaxEnt production run...");

    auto job = mdrun()
        .multidir(rep_dirs)
        .plumed(plumed_file)
        .run_multi(std::move(prod_runner));
    job->wait();

    //##################################//
    //###      POST-PROCESSING       ###//
    //##################################//
    // Strip solvent and re-center each replica trajectory.
    console::print_text("Post-processing replica trajectories...");
    for (int i = 0; i < N; ++i) {
        const auto& rep_dir = rep_dirs[i];
        GROFile prodgro(rep_dir + "prod.gro");
        XTCFile prodxtc(rep_dir + "prod.xtc");
        TPRFile prodtpr(rep_dir + "prod.tpr");

        auto [ndx_protein] = md::select(prodgro)
            .group("Protein")
            .output(rep_dir + "protein.ndx")
        .run();

        trjconv(prodxtc)
            .group("Protein")   // echo "Protein 0 |" for center + output groups
            .center()
            .pbc("mol")
            .runfile(prodtpr)
            .index(ndx_protein)
            .output(rep_dir + "protein.xtc")
        .run();

        editconf(prodgro)
            .index(ndx_protein)
            .output(rep_dir + "protein.gro")
        .run();

        ndx_protein.remove();
    }

    console::print_info("Done. Results in " + settings::general::output
        + "\n  Trajectories : prod/rep{i}/protein.xtc"
        + "\n  PLUMED output: prod/rep{i}/COLVAR");
    return 0;
}
