#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/frameanalysis.h>

#include <CLI/CLI.hpp>

using namespace gmx;

int main(int argc, char const *argv[]) {
    std::string s_pdb;
    CLI::App app{"Perform an analysis of the SAXS-dependency on the number of frames."};
    app.add_option("input", s_pdb, "PDB structure file.")->required();
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
        .output = "output/frame_analysis/" + std::filesystem::path(s_pdb).stem().string() + "/",
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute(),
        .setupsim = location::lucy,
        .mainsim = location::smaug,
        .bufmdp = std::make_shared<FrameAnalysisMDPCreatorSol>(),
        .molmdp = std::make_shared<FrameAnalysisMDPCreatorMol>(),
    };
    gmx::gmx::set_cmdlog(sele.output + "cmd.log");
    gmx::gmx::set_outputlog(sele.output + "output.log");
    PDBFile pdb(s_pdb);

    // prepare sims
    MoleculeOptions mo(sele, pdb);
    auto molecule = simulate_molecule(mo);

    gmx::SimulateBufferOutput buffer;
    BufferOptions bo(sele, molecule.gro);
    buffer = simulate_buffer(bo);

    // prepare saxs
    GMXOptions saxs_sele = sele;
    saxs_sele.jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute();
    
    SAXSOptions so(saxs_sele, std::move(molecule), std::move(buffer), pdb);
    auto saxs = frameanalysis(so);

    for (unsigned int i = 0; i < saxs.size(); i++) {saxs[i].job->submit();}
    for (unsigned int i = saxs.size(); i > 0; i--) { // reverse wait order since last job is the shortest one
        saxs[i].job->wait();
        auto res = saxs[i].job->result();
        res.xvg.rename("saxs_" + std::to_string(i) + ".xvg").copy(sele.output + "output");
    }

    return 0;
} 