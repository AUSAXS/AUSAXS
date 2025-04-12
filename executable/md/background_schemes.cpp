#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/frameanalysis.h>
#include <io/ExistingFile.h>

#include <CLI/CLI.hpp>

using namespace ausaxs;
using namespace ausaxs::md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb;
    CLI::App app{"Calculate two scattering curves using both background subtraction schemes for a single MD simulation."};
    app.add_option("input", s_pdb, "PDB structure file.")->required();
    CLI11_PARSE(app, argc, argv);

    // test executable
    if (!gmx().valid_executable()) {
        throw except::io_error("Gromacs executable not found. Please install Gromacs and add it to your PATH.");
    }
    
    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = s_pdb.stem(),
        .output = {"output/frame_analysis/" + s_pdb.stem() + "/"},
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
        .setupsim = RunLocation::lusi,
        .mainsim = RunLocation::smaug,
        .bufmdp = std::make_shared<FrameAnalysisMDPCreatorSol>(),
        .molmdp = std::make_shared<FrameAnalysisMDPCreatorMol>(),
    };

    if (io::Folder tmp("temp/md"); !tmp.exists()) {tmp.create();}
    gmx::gmx::set_logfile(sele.output.str() + "output.log", sele.output.str() + "cmd.log");
    PDBFile pdb(s_pdb);

    // prepare sims
    MoleculeOptions mo(sele, pdb);
    auto molecule = simulate_molecule(mo);

    SimulateBufferOutput buffer;
    BufferOptions bo(sele, molecule.gro);
    buffer = simulate_buffer(bo);

    // prepare saxs
    GMXOptions saxs_sele = sele;
    saxs_sele.jobscript = SHFile("temp/scripts/jobscript_slurm_swaxs.sh").absolute_path();
    
    SAXSOptions so(saxs_sele, std::move(molecule), std::move(buffer), pdb);
    auto saxs = frameanalysis(so);

    for (unsigned int i = 0; i < saxs.size(); i++) {saxs[i].job->submit();}
    for (unsigned int i = saxs.size(); i > 0; i--) { // reverse wait order since last job is the shortest one
        saxs[i].job->wait();
        auto res = saxs[i].job->result();
        res.xvg.rename("saxs_" + std::to_string(i) + ".xvg").copy({sele.output.str() + "output"});
    }

    return 0;
}