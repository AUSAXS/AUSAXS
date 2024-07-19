#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/saxs.h>
#include <settings/MDSettings.h>
#include <settings/SettingsIO.h>
#include <io/ExistingFile.h>
#include <constants/Version.h>

#include <CLI/CLI.hpp>

using namespace md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb, s_settings;
    CLI::App app{"MD simulation pipeline."};
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
        if (settings::discover(".")) {
            CLI11_PARSE(app, argc, argv);
        }
    }
    if (!gmx().valid_executable()) {
        throw except::io_error("GROMACS path \"" + settings::md::gmx_path + "\"");
    }

    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = s_pdb.stem(),
        .output = "output/md/" + s_pdb.stem() + "/",
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute(),
        .setupsim = location::local,
        .mainsim = location::local,
        .bufmdp = std::make_shared<PRMDPCreatorSol>(),
        .molmdp = std::make_shared<PRMDPCreatorMol>(),
    };

    // gmx::gmx::set_cmdlog(sele.output + "cmd.log");
    // gmx::gmx::set_outputlog(sele.output + "output.log");
    PDBFile pdb(s_pdb);

    // prepare sims
    MoleculeOptions mo(sele, pdb);
    auto molecule = simulate_molecule(mo);

    SimulateBufferOutput buffer;
    if (settings::md::buffer_path.empty()) {
        BufferOptions bo(sele, molecule.gro);
        buffer = simulate_buffer(bo);
    } else {
        console::print_info("\nPreparing buffer simulation");
        // find production gro
        for (auto& p : std::filesystem::directory_iterator(settings::md::buffer_path + "/prod")) {
            if (p.path().extension() == ".gro") {
                buffer.gro = GROFile(p.path());
                std::cout << "\tFound production gro file: " << buffer.gro.path << std::endl;
            }
        }

        // find topology file
        for (auto& p : std::filesystem::directory_iterator(settings::md::buffer_path + "/setup")) {
            if (p.path().extension() == ".top") {
                buffer.top = TOPFile(p.path());
                std::cout << "\tFound topology file: " << buffer.top.path << std::endl;
            }
        }

        // create dummy job
        buffer.job = std::make_unique<NoExecution<MDRunResult>>(settings::md::buffer_path + "/prod/");

        // check that all files are found
        if (buffer.gro.empty() || buffer.top.empty()) {
            throw except::io_error("Could not find all files in supplied buffer folder \"" + settings::md::buffer_path + "\".");
        }
    }

    // prepare saxs
    GMXOptions saxs_sele = sele;
    saxs_sele.jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute();
    
    SAXSOptions so(saxs_sele, std::move(molecule), std::move(buffer), pdb);
    auto saxs = simulate_saxs(so);
    saxs.job->submit();

    return 0;
} 