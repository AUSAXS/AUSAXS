/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/simulate/SimulateMolecule.h>
#include <utility/Console.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;
using namespace ausaxs::md;

SimulateMoleculeOutput md::simulate_molecule(SimulateMoleculeOptions&& options) {
    console::print_info("\nPreparing simulation for " + options.pdbfile.filename());
    console::indent();

    //##################################//
    //###           GLOBALS          ###//
    //##################################//
    io::Folder mdp_folder = settings::general::output + "mdp/";
    io::Folder setup_path = settings::general::output + "protein/setup/";
    io::Folder em_path    = settings::general::output + "protein/em/";
    io::Folder eq_path    = settings::general::output + "protein/eq/";
    io::Folder prod_path  = settings::general::output + "protein/prod/";
    mdp_folder.create(); setup_path.create(); em_path.create(); eq_path.create(); prod_path.create();

    // gmx::gmx::set_outputlog(output + "gmx.log");
    //##################################//
    //###           SETUP            ###//
    //##################################//
    GROFile solv_ion(setup_path + "solv_ion.gro");
    TOPFile top(setup_path + "topol.top");
    auto ff = option::IForcefield::construct(options.system.forcefield);
    auto wm = option::IWaterModel::construct(options.system.watermodel);
    if (!solv_ion.exists() || !top.exists()) {
        GROFile conf(setup_path + "conf.gro");
        if (!conf.exists()) {
            console::print_text("Converting structure to GROMACS format...");

            // prepare the pdb file for gromacs
            auto cmd = pdb2gmx(options.pdbfile)
                .output(setup_path)
                .ignore_hydrogens()
                .water_model(wm.get())
                .forcefield(ff.get());
            if (options.system.custom_options.contains("vsites") && options.system.custom_options.at("vsites") == "true") {cmd.virtual_sites();}
            auto[conf, top, posre] = cmd.run();
        } else {
            console::print_text("Reusing previously generated GROMACS structure.");
        }

        // create a box around the protein
        console::print_text("Generating unit cell...");
        auto[uc] = editconf(conf)
            .output(setup_path + "uc.gro")
            .box_type(options.system.boxtype)
            .extend(options.system.custom_options.contains("editconf extend") ? std::stod(options.system.custom_options.at("editconf extend")) : 1.5)
        .run();

        // add water to the box
        console::print_text("Solvating unit cell...");
        auto[solv] = solvate(uc)
            .output(setup_path + "solv.gro")
            .solvent(ff.get(), wm.get())
            .radius(0.2)
            .topology(top)
        .run();

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
            .cation(options.system.cation)
            .anion(options.system.anion)
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

        // prepare energy minimization sim
        MDPFile mdp = mdp::templates::minimize::base().write(mdp_folder + "emmol.mdp");
        auto[emtpr] = grompp(mdp, top, solv_ion)
            .output(em_path + "em.tpr")
            .restraints(solv_ion)
        .run();

        // run energy minimization
        mdrun(emtpr)
            .output(em_path, "/em")
        .run(std::move(options.minimize_runner))->wait();
    } else {
        console::print_text("Reusing previously generated energy minimization.");
    }

    //##################################//
    //###       EQUILIBRATION        ###//
    //##################################//
    GROFile eqgro(eq_path + "eq.gro");
    NDXFile index(eq_path + "index.ndx");
    if (!eqgro.exists()) {
        console::print_text("Running thermalization...");

        // make a basic index file without any special groups
        auto[_] = make_ndx(emgro)
            .output(index)
        .run();

        // thermalization requires the group "Water_and_Ions". 
        // If it is not present, it is because we have no ions, so it is just the water group
        if (!index.contains("Water_and_Ions")) {
            auto[tmp] = select(emgro)
                .output("tmp.ndx")
                .define("\"Water_and_Ions\" group \"Water\"")
            .run();

            index.append_file(tmp);
            tmp.remove();
        }

        // if (!index.contains("Protein")) {
        //     std::string group_def;
        //     if (index.contains("DNA")) {group_def += " group \"DNA\"";}
        //     if (index.contains("RNA")) {
        //         if (!group_def.empty()) {group_def += " or ";}
        //         group_def += "group \"RNA\"";
        //     }
        //     if (group_def.empty()) {throw std::runtime_error("molecule: No DNA or RNA found in the index file");}

        //     // add a protein group to the index file
        //     auto[tmp] = select(emgro)
        //         .output("tmp.ndx")
        //         .define("\"Protein\"" + group_def)
        //     .run();

        //     index.append_file(tmp);
        //     tmp.remove();
        // }

        // prepare equilibration sim
        MDPFile mdp = mdp::templates::equilibrate::mol().write(mdp_folder + "eqmol.mdp");
        auto[eqtpr] = grompp(mdp, top, emgro)
            .output(eq_path + "eq.tpr")
            .index(index)
            .restraints(emgro)
            .warnings(1)
        .run();

        // run equilibration
        mdrun(eqtpr)
            .output(eq_path, "/eq")
        .run(std::move(options.equilibrate_runner))->wait();
    } else {
        console::print_text("Reusing previously generated thermalization.");
    }

    //##################################//
    //###       PRODUCTION           ###//
    //##################################//
    GROFile prodgro(prod_path + "prod.gro");
    if (!prodgro.exists()) {
        if (options.system.custom_options.contains("backbone_restraints") && options.system.custom_options.at("backbone_restraints") == "true") {
            console::print_text("Generating backbone restraints...");
            auto[backbone] = genrestr(eq_path + "eq.gro")
                .output(setup_path + "backbone.itp")
                .index(index)
                .force(2000, 2000, 2000)
            .run();
            auto backbone_chains = backbone.split_restraints(top.get_includes());
            top.include_new_type(backbone_chains, "POSRESBACKBONE");
        } else {
            console::print_text("No backbone restraints will be used.");
        }
        // top.include(top.relative_path(backbone), "POSRESBACKBONE", "chain topol");

        // prepare production sim
        console::print_text("Running production...");
        auto cmd = grompp(options.mdp, top, eqgro)
            .output(prod_path + "prod.tpr")
            .index(index)
            .warnings(1);
        if (options.system.custom_options.contains("backbone_restraints") && options.system.custom_options.at("backbone_restraints") == "true") {cmd.restraints(eqgro);}
        auto[prodtpr] = cmd.run();

        // run production
        auto job = mdrun(prodtpr)
            .output(prod_path, "/prod")
        .run(std::move(options.production_runner));

        console::unindent();
        return {std::move(job), top, solv_ion};
    } else {
        console::print_text("Reusing previous position-restrained simulation.");
    }
    console::unindent();
    auto job = std::make_unique<NoExecutor<MDRunResult>>(prod_path);
    return {std::move(job), top, solv_ion};
}