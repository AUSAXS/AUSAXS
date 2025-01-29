/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/simulate/molecule.h>
#include <utility/Console.h>

using namespace ausaxs;
using namespace ausaxs::md;

md::simulate::Molecule::Molecule(MoleculeOptions& options) : options(options) {}

md::SimulateMoleculeOutput md::simulate::Molecule::simulate() {return SimulateMoleculeOutput();}

std::tuple<GROFile, TOPFile> md::simulate::Molecule::setup() {return std::make_tuple(GROFile(), TOPFile());}

GROFile md::simulate::Molecule::minimize() {return GROFile();}

std::tuple<GROFile, NDXFile> md::simulate::Molecule::thermalize() {return std::make_tuple(GROFile(), NDXFile());}

SimulateMoleculeOutput md::simulate_molecule(MoleculeOptions& options) {
    console::print_info("\nPreparing simulation for " + options.pdbfile.filename());
    console::indent();

    //##################################//
    //###           GLOBALS          ###//
    //##################################//
    io::Folder mdp_folder = options.output.str() + "mdp/";
    io::Folder setup_path = options.output.str() + "protein/setup/";
    io::Folder em_path = options.output.str() + "protein/em/";
    io::Folder eq_path = options.output.str() + "protein/eq/";
    io::Folder prod_path = options.output.str() + "protein/prod/";
    mdp_folder.create(); setup_path.create(); em_path.create(); eq_path.create(); prod_path.create();

    // gmx::gmx::set_outputlog(output + "gmx.log");
    //##################################//
    //###           SETUP            ###//
    //##################################//
    GROFile solv_ion(setup_path.str() + "solv_ion.gro");
    TOPFile top(setup_path.str() + "topol.top");
    if (!solv_ion.exists() || !top.exists()) {
        GROFile conf(setup_path.str() + "conf.gro");
        if (!conf.exists()) {
            console::print_text("Converting structure to GROMACS format...");

            // prepare the pdb file for gromacs
            auto[conf, top, posre] = pdb2gmx(options.pdbfile)
                .output(setup_path)
                .ignore_hydrogens()
                // .virtual_sites()
                .water_model(options.watermodel)
                .forcefield(options.forcefield)
            .run();
        } else {
            console::print_text("Reusing previously generated GROMACS structure.");
        }

        // create a box around the protein
        console::print_text("Generating unit cell...");
        auto[uc] = editconf(conf)
            .output(setup_path.str() + "uc.gro")
            .box_type(options.boxtype)
            .extend(1.5)
        .run();

        // add water to the box
        console::print_text("Solvating unit cell...");
        auto[solv] = solvate(uc)
            .output(setup_path.str() + "solv.gro")
            .solvent(options.forcefield, options.watermodel)
            .radius(0.2)
            .topology(top)
        .run();

        // generate an empty tpr file 
        MDPFile mdp = MDPFile(setup_path.str() + "empty.mdp"); mdp.create();
        auto[ions] = grompp(mdp, top, solv)
            .output(setup_path.str() + "ions.tpr")
            .warnings(1)
        .run();
        mdp.remove();

        // add ions to the box
        console::print_text("Adding ions to unit cell...");
        genion(ions)
            .output(solv_ion)
            .topology(top)
            .neutralize()
            .cation(options.cation)
            .anion(options.anion)
            .ion_group("SOL")
        .run();
    } else {
        console::print_text("Reusing previously generated system setup.");
    }

    //##################################//
    //###     ENERGY MINIMIZATION    ###//
    //##################################//
    GROFile emgro(em_path.str() + "em.gro");
    if (!emgro.exists()) {
        console::print_text("Running energy minimization...");

        // prepare energy minimization sim
        MDPFile mdp = EMMDPCreator().write(mdp_folder.str() + "emmol.mdp");
        auto[emtpr] = grompp(mdp, top, solv_ion)
            .output(em_path.str() + "em.tpr")
            .restraints(solv_ion)
        .run();

        // run energy minimization
        mdrun(emtpr)
            .output(em_path, "em")
            .jobname(options.name + "_mol")
        .run(options.setupsim, options.jobscript)->submit();
    } else {
        console::print_text("Reusing previously generated energy minimization.");
    }

    //##################################//
    //###       EQUILIBRATION        ###//
    //##################################//
    GROFile eqgro(eq_path.str() + "eq.gro");
    NDXFile index(eq_path.str() + "index.ndx");
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
        MDPFile mdp = EQMDPCreatorMol().write(mdp_folder.str() + "eqmol.mdp");
        auto[eqtpr] = grompp(mdp, top, emgro)
            .output(eq_path.str() + "eq.tpr")
            .index(index)
            .restraints(emgro)
            .warnings(1)
        .run();

        // run equilibration
        mdrun(eqtpr)
            .output(eq_path, "eq")
            .jobname(options.name + "_mol")
        .run(options.setupsim, options.jobscript)->submit();
    } else {
        console::print_text("Reusing previously generated thermalization.");
    }

    //##################################//
    //###       PRODUCTION           ###//
    //##################################//
    GROFile prodgro(prod_path.str() + "prod.gro");
    if (!prodgro.exists()) {
        console::print_text("Generating backbone restraints...");
        auto[backbone] = genrestr(eq_path.str() + "eq.gro")
            .output(setup_path.str() + "backbone.itp")
            .index(index)
            .force(2000, 2000, 2000)
        .run();
        auto backbone_chains = backbone.split_restraints(top.get_includes());
        top.include_new_type(backbone_chains, "POSRESBACKBONE");
        // top.include(top.relative_path(backbone), "POSRESBACKBONE", "chain topol");

        // prepare production sim
        console::print_text("Running production...");
        MDPFile mdp = options.molmdp->write(mdp_folder.str() + "prmol.mdp");
        auto[prodtpr] = grompp(mdp, top, eqgro)
            .output(prod_path.str() + "prod.tpr")
            .index(index)
            .restraints(eqgro)
            .warnings(1)
        .run();

        // run production
        auto job = mdrun(prodtpr)
            .output(prod_path, "prod")
            .jobname(options.name + "_mol")
        .run(options.mainsim, options.jobscript);

        console::unindent();
        return {std::move(job), top, solv_ion};
    } else {
        console::print_text("Reusing previous position-restrained simulation.");
    }
    console::unindent();
    auto job = std::make_unique<NoExecution<MDRunResult>>(prod_path);
    return {std::move(job), top, solv_ion};
}