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
#include <dataset/SimpleDataset.h>
#include <plots/All.h>

#include <CLI/CLI.hpp>

using namespace ausaxs;
using namespace ausaxs::md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb, s_settings;
    CLI::App app{"Comparison of the two background subtraction schemes."};
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
    settings::general::output += "md/" + s_pdb.stem() + "/";

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(settings::general::output + "output.log", settings::general::output + "cmd.log");
    PDBFile pdb(s_pdb);

    SystemSettings ss {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,
    };

    auto molecule = simulate_molecule({
        .system = ss,
        .jobname = pdb.stem() + "_mol",
        .pdbfile = pdb,
        .mdp = mdp::templates::production::mol().write(settings::general::output + "mdp/mol.mdp"),
        .setup_runner = RunLocation::lusi,
        .main_runner = RunLocation::smaug,
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
    });

    SimulateBufferOutput buffer;
    console::print_info("\nPreparing buffer simulation");
    if (settings::md::buffer_path.empty()) {
        buffer = simulate_buffer({
            .system = ss,
            .jobname = pdb.stem() + "_buf",
            .refgro = molecule.gro,
            .mdp = mdp::templates::production::solv().write(settings::general::output + "mdp/buf.mdp"),
            .setup_runner = RunLocation::smaug,
            .main_runner = RunLocation::smaug,
            .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
        });
    } else {
        // find production gro
        io::Folder prod_folder(settings::md::buffer_path + "/prod");
        if (!prod_folder.exists()) {
            throw except::io_error("Could not find nested production folder (\"prod\") in supplied buffer folder \"" + settings::md::buffer_path + "\".");
        }
        for (auto& p : prod_folder.files()) {
            if (p.extension() == ".gro") {
                buffer.gro = GROFile(p.path());
                console::print_text("\tFound production gro file: " + buffer.gro.path());
            }
        }
        if (!buffer.gro.exists()) {
            throw except::io_error("Could not find production gro file in supplied buffer folder \"" + settings::md::buffer_path + "/prod/\".");
        }

        // find topology file
        io::Folder setup_folder(settings::md::buffer_path + "/setup");
        for (auto& p : setup_folder.files()) {
            if (p.extension() == ".top") {
                buffer.top = TOPFile(p.path());
                console::print_text("\tFound topology file: " + buffer.top.path());
            }
        }
        if (!buffer.top.exists()) {
            throw except::io_error("Could not find topology file in supplied buffer folder \"" + settings::md::buffer_path + "/setup/\".");
        }

        // create dummy job
        buffer.job = std::make_unique<NoExecution<MDRunResult>>(settings::md::buffer_path + "/prod/");
    }

    SimulateSAXSOutput saxs1, saxs2;
    {   // scheme 1: uncorrected
        auto mol_mdp = mdp::templates::saxs::mol();
        auto buf_mdp = mdp::templates::saxs::solv();
        mol_mdp.add(MDPOptions::waxs_correct_buffer = "no");
        buf_mdp.add(MDPOptions::waxs_correct_buffer = "no");

        saxs1 = simulate_saxs({
            .jobname = "" + pdb.stem() + "_saxs1",
            .pdbfile = pdb,
            .molecule = molecule,
            .buffer = buffer,
            .runner = RunLocation::smaug,
            .jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute_path(),
            .output = settings::general::output + "saxs_uncorrected/",
            .mol_mdp = mol_mdp,
            .buf_mdp = buf_mdp,
        });
        saxs1.job->submit();
    }
    {   // scheme 2: corrected
        auto mol_mdp = mdp::templates::saxs::mol();
        auto buf_mdp = mdp::templates::saxs::solv();
        mol_mdp.add(MDPOptions::waxs_correct_buffer = "yes");
        buf_mdp.add(MDPOptions::waxs_correct_buffer = "yes");

        saxs2 = simulate_saxs({
            .jobname = "" + pdb.stem() + "_saxs2",
            .pdbfile = pdb,
            .molecule = molecule,
            .buffer = buffer,
            .runner = RunLocation::smaug,
            .jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute_path(),
            .output = settings::general::output + "saxs_corrected/",
            .mol_mdp = mol_mdp,
            .buf_mdp = buf_mdp,
        });
        saxs2.job->submit();
    }

    SimpleDataset data1(saxs1.job->result().xvg);
    SimpleDataset data2(saxs2.job->result().xvg);
    plots::PlotDataset()
        .plot(data1, {style::draw::line, {{"legend", "uncorrected"}, {"logx", true}, {"logy", true}, {"color", style::color::orange}}})
        .plot(data2, {style::draw::line, {{"legend", "corrected"}, {"color", style::color::blue}}})
    .save(settings::general::output + "out/saxs_comparison.png");
    saxs1.job->result().xvg.copy(settings::general::output + "out/saxs_uncorrected.xvg");
    saxs2.job->result().xvg.copy(settings::general::output + "out/saxs_corrected.xvg");
    return 0;
}