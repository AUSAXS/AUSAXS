#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/saxs.h>
#include <md/gmx/Settings.h>

#include <CLI/CLI.hpp>

using namespace gmx;

int main(int argc, char const *argv[]) {
    settings::discover(".");

    std::string s_pdb;
    CLI::App app{"MD simulation pipeline."};
    app.add_option("input", s_pdb, "PDB structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--buffer", setting::buffer_path.value, "Pre-simulated buffer. Set this to the base output folder.");
    CLI11_PARSE(app, argc, argv);

    // test executable
    if (!gmx::gmx().test_executable()) {
        throw except::io_error("Gromacs executable not found. Please install Gromacs and add it to your PATH.");
    }

    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = std::filesystem::path(s_pdb).stem().string(),
        .output = "output/" + std::filesystem::path(s_pdb).stem().string() + "/",
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

    gmx::SimulateBufferOutput buffer;
    if (setting::buffer_path.get().empty()) {
        BufferOptions bo(sele, molecule.gro);
        buffer = simulate_buffer(bo);
    } else {
        utility::print_info("\nPreparing buffer simulation");
        // find production gro
        for (auto& p : std::filesystem::directory_iterator(setting::buffer_path.get() + "/prod")) {
            if (p.path().extension() == ".gro") {
                buffer.gro = GROFile(p.path());
                std::cout << "\tFound production gro file: " << buffer.gro.path << std::endl;
            }
        }

        // find topology file
        for (auto& p : std::filesystem::directory_iterator(setting::buffer_path.get() + "/setup")) {
            if (p.path().extension() == ".top") {
                buffer.top = TOPFile(p.path());
                std::cout << "\tFound topology file: " << buffer.top.path << std::endl;
            }
        }

        // create dummy job
        buffer.job = std::make_unique<NoExecution<MDRunResult>>(setting::buffer_path.get() + "/prod/");

        // check that all files are found
        if (buffer.gro.empty() || buffer.top.empty()) {
            throw except::io_error("Could not find all files in supplied buffer folder \"" + setting::buffer_path.get() + "\".");
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