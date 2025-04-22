/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/simulate/SimulateSAXS.h>
#include <md/programs/all.h>
#include <md/programs/saxsmdrun.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/utility/Exceptions.h>
#include <md/utility/Protein.h>
#include <md/utility/files/MDPCreator.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>
#include <io/Folder.h>

#include <cmath>

using namespace ausaxs;

md::SimulateSAXSOutput md::simulate_saxs(md::SimulateSAXSOptions&& options) {
    if (!options.molecule.top.exists()) {throw except::io_error("simulate_saxs: The topology file does not exist.");}

    // get data from sims
    if (options.molecule.job == nullptr) {throw except::unexpected("gmx::simulate: No molecule simulation was created.");}
    if (options.buffer.job == nullptr) {throw except::unexpected("gmx::simulate: No buffer was created.");}
    options.molecule.job->wait();
    options.buffer.job->wait();
    auto[molgro, moledr, molxtc] = options.molecule.job->result();
    auto[bufgro, bufedr, bufxtc] = options.buffer.job->result();

    console::print_info("\nPreparing calculation of SAXS profile");
    console::indent();
    //##################################//
    //###           GLOBALS          ###//
    //##################################//
    io::Folder protein_path = options.output + "protein/";
    io::Folder buffer_path  = options.output + "buffer/";
    io::Folder mdp_folder   = options.output + "mdp/";
    io::Folder prod_folder  = options.output + "prod/";
    mdp_folder.create(); prod_folder.create();

    TOPFile moltop(protein_path + "topol.top");
    TOPFile buftop(buffer_path + "topol.top");
    if (!moltop.exists() || !buftop.exists()) {
        options.molecule.top.copy(moltop.directory());
        options.buffer.top.copy(buftop.directory());
    } else {
        console::print_text("Reusing topology files for molecule and buffer.");
    }

    std::cout << "\tChecking if trajectories are compatible..." << std::flush;
    {
        auto[molframes, molduration] = check(molxtc).run();
        auto[bufframes, bufduration] = check(bufxtc).run();
        if (molframes != bufframes) {throw except::unexpected("simulate_saxs: The number of frames in the buffer and molecule trajectories are not equal.");}
        if (molduration != bufduration) {throw except::unexpected("simulate_saxs: The duration of the buffer and molecule trajectories are not equal.");}
    }
    std::cout << " OK" << std::endl;

    //##################################//
    //###        INDEX FILES         ###//
    //##################################//
    // we want to always have the group Water_and_ions, but it is only generated automatically if there are actual ions present in the system. 
    // since this is not guaranteed, we have to generate it manually
    NDXFile molindex(protein_path + "index.ndx");
    if (!molindex.exists()) {
        console::print_text("Creating index file for molecule...");
        make_ndx(molgro)
            .output(molindex)
        .run();

        std::string ion_placeholder = molindex.contains("Ion") ? " or group \"Ion\"" : "";
        if (!molindex.contains("RealWater_and_Ions")) {
            // if there are no Ion group, adding "or group \"Ion\" will cause an error
            auto[tmp] = select(options.molecule.gro)
                .output("temp/md/tmp1.ndx")
                .define("\"RealWater_and_Ions\" name \"OW\" or name \"HW1\" or name \"HW2\" or name \"HW3\"" + ion_placeholder)
            .run();

            molindex.append_file(tmp);
            tmp.remove();
        }

        if (!molindex.contains("Water_and_Ions")) {
            auto[tmp] = select(options.molecule.gro)
                .output("temp/md/tmp2.ndx")
                .define("\"Water_and_Ions\" name \"OW\" or name \"HW1\" or name \"HW2\" or name \"HW3\"" + ion_placeholder)
            .run();

            molindex.append_file(tmp);
            tmp.remove();
        }

        if (!molindex.contains("Protein_and_Ions")) {
            auto[tmp] = select(options.molecule.gro)
                .output("temp/md/tmp3.ndx")
                .define("\"Protein_and_Ions\" group \"Protein\"" + ion_placeholder)
            .run();

            molindex.append_file(tmp);
            tmp.remove();
        }
    } else {
        console::print_text("Reusing index file for molecule.");
    }

    NDXFile bufindex(buffer_path + "index.ndx");
    if (!bufindex.exists()) {
        console::print_text("Creating index file for buffer...");
        make_ndx(bufgro)
            .output(bufindex)
        .run();
        if (!bufindex.contains("RealWater")) {
            auto[tmp] = select(options.buffer.gro)
                .output("tmp.ndx")
                .define("\"RealWater\" name \"OW\" or name \"HW1\" or name \"HW2\" or name \"HW3\"")
            .run();

            bufindex.append_file(tmp);
            tmp.remove();
        }
    } else {
        console::print_text("Reusing index file for buffer.");
    }

    //##################################//
    //###         ENVELOPE           ###//
    //##################################//
    GROFile envgro(protein_path + "envelope-ref.gro");
    PYFile envpy(protein_path + "envelope.py");
    DATFile envdat(protein_path + "envelope.dat");
    MDPFile molmdp(mdp_folder + "rerun_mol.mdp");

    if (!envgro.exists() || !envpy.exists() || !envdat.exists() || !molmdp.exists()) {
        console::print_text("Preparing SAXS rerun...");
        console::indent();
        MDPFile dummymdp = MDPFile(options.output + "empty.mdp"); dummymdp.create();
        auto[dummytpr] = grompp(dummymdp, moltop, molgro)
            .output(options.output + "saxs.tpr")
            .warnings(1)
        .run();
        dummymdp.remove();

        console::print_text("Generating scattering parameters...");
        auto[itps] = genscatt(dummytpr, molindex)
            .output(protein_path + "scatt.itp")
            .group("Protein_and_Ions")
        .run();

        // remove NA and CL scattering parameters
        for (unsigned int i = 0; i < itps.size(); i++) {
            auto itp = itps[i].stem().substr(0, 8);
            if (itp == "scatt_NA" || itp == "scatt_CL") {
                itps[i].remove();
                itps.erase(itps.begin() + i);
            }
        }

        moltop.include_new_type(itps);

        // dump the first 50 frames
        console::print_text("Dumping first 50 frames...");
        auto[traj] = trjconv(molxtc)
            .output(protein_path + "protein.xtc")
            .startframe(50)
        .run();

        // center the trajectories
        console::print_text("Centering trajectories...");
        auto[cluster] = trjconv(traj)
            .output(protein_path + "cluster.xtc")
            .group("Protein")
            .pbc("cluster")
            .ur("tric")
            .index(molindex)
            .runfile(dummytpr)
        .run();

        auto[centered] = trjconv(cluster)
            .output(protein_path + "centered.xtc")
            .group("Protein")
            .center()
            .boxcenter("tric")
            .pbc("mol")
            .ur("tric")
            .index(molindex)
            .runfile(dummytpr)
        .run();

        // generate the envelope
        console::print_text("Generating envelope...");
        auto[goodpbc, envgro, envpy, envdat] = genenv(centered, molindex)
            .output(protein_path)
            .structure(molgro)
            .distance(0.5)
            .group("Protein")
        .run();
        dummytpr.remove();

        // generate molecule mdp file (depends on envelope output)
        {
            auto _mdp = options.mol_mdp.value_or(mdp::templates::saxs::mol());

            Protein protein(options.pdbfile);
            double qmax = std::stod(_mdp.get(MDPOptions::waxs_endq))/10; // convert to nm^-1
            double dmax = protein.maximum_distance();

            _mdp.add(MDPOptions::waxs_pbcatom = goodpbc);
            _mdp.add(MDPOptions::waxs_nsphere = int(0.2*std::pow((dmax*qmax), 2)));
            molmdp = _mdp.write(molmdp);
        }
        console::unindent();
    } else {
        console::print_text("Reusing previously generated envelope.");
    }

    MDPFile bufmdp(mdp_folder + "rerun_buf.mdp");
    if (!bufmdp.exists()) {
        options.buf_mdp.value_or(mdp::templates::saxs::solv()).write(bufmdp);
    }

    GROFile gro(prod_folder + "prod.gro");
    if (!gro.exists()) {
        console::print_text("Rerunning production simulation...");
        auto[moltpr] = grompp(molmdp, moltop, molgro)
            .output(prod_folder + "rerun_mol.tpr")
            .index(molindex)
            .warnings(1)
        .run();

        auto[buftpr] = grompp(bufmdp, buftop, bufgro)
            .output(prod_folder + "rerun_buf.tpr")
            .index(bufindex)
            .warnings(2)
        .run();

        auto job = saxsmdrun(moltpr, buftpr)
            .output(prod_folder, "prod")
            .rerun(molxtc, bufxtc)
            .env_var("GMX_WAXS_FIT_REFFILE", envgro.absolute_path())
            .env_var("GMX_ENVELOPE_FILE", envdat.absolute_path())
            .jobname(options.jobname)
        .run(options.runner, options.jobscript);

        console::unindent();
        return SimulateSAXSOutput{std::move(job)};
    } else {
        console::print_text("Reusing previously generated final simulation.");
    }
    console::unindent();
    auto job = std::make_unique<NoExecution<SAXSRunResult>>(prod_folder);
    return SimulateSAXSOutput{std::move(job)};
}