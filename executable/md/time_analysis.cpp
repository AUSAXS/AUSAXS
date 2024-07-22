#include <md/programs/all.h>
#include <md/simulate/buffer.h>
#include <md/simulate/molecule.h>
#include <md/simulate/timeanalysis.h>
#include <io/ExistingFile.h>

#include <CLI/CLI.hpp>

using namespace md;

int main(int argc, char const *argv[]) {
    io::ExistingFile s_pdb;
    CLI::App app{"Perform an analysis of the SAXS-dependency on the simulation time."};
    app.add_option("input", s_pdb, "PDB structure file.")->required();
    CLI11_PARSE(app, argc, argv);

    // test executable
    if (!gmx().valid_executable()) {
        throw except::io_error("Gromacs executable not found. Please install Gromacs and add it to your PATH.");
    }

    GMXOptions sele {
        .forcefield = option::Forcefield::AMBER99SB_ILDN,
        .watermodel = option::WaterModel::TIP4P2005,
        .boxtype = option::BoxType::DODECAHEDRON,
        .cation = option::Cation::NA,
        .anion = option::Anion::CL,

        .name = s_pdb.stem(),
        .output = std::string_view("output/time_analysis/" + s_pdb.stem() + "/"),
        .jobscript = SHFile("scripts/jobscript_slurm_standard.sh").absolute_path(),
        .setupsim = location::lucy,
        .mainsim = location::smaug,
        .bufmdp = std::make_shared<TimeAnalysisMDPCreatorSol>(),
        .molmdp = std::make_shared<TimeAnalysisMDPCreatorMol>(),
    };
    gmx::gmx::set_cmdlog(sele.output + "cmd.log");
    gmx::gmx::set_outputlog(sele.output + "output.log");
    PDBFile pdb(s_pdb);

    // prepare sims
    MoleculeOptions mo(sele, pdb);
    auto molecule = simulate_molecule(mo);

    SimulateBufferOutput buffer;
    BufferOptions bo(sele, molecule.gro);
    buffer = simulate_buffer(bo);

    // prepare saxs
    GMXOptions saxs_sele = sele;
    saxs_sele.jobscript = SHFile("scripts/jobscript_slurm_swaxs.sh").absolute_path();
    
    SAXSOptions so(saxs_sele, std::move(molecule), std::move(buffer), pdb);
    auto saxs = timeanalysis(so, 500);

    for (unsigned int i = 0; i < saxs.size(); i++) {saxs[i].job->submit();}
    for (unsigned int i = saxs.size()-1; i >= 0; i--) { // reverse wait order since the last job is the shortest run
        saxs[i].job->wait();
        auto res = saxs[i].job->result();
        res.xvg.rename(std::to_string(i+1) + ".xvg").copy({sele.output + "output"});
    }
    return 0;
} 