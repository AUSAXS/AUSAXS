/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/simulate/molecule.h>
#include <utility/Console.h>

using namespace md;

md::simulate::Molecule::Molecule(MoleculeOptions& options) : options(options) {}

md::SimulateMoleculeOutput md::simulate::Molecule::simulate() {return SimulateMoleculeOutput();}

std::tuple<GROFile, TOPFile> md::simulate::Molecule::setup() {return std::make_tuple(GROFile(), TOPFile());}

GROFile md::simulate::Molecule::minimize() {return GROFile();}

std::tuple<GROFile, NDXFile> md::simulate::Molecule::thermalize() {return std::make_tuple(GROFile(), NDXFile());}

SimulateMoleculeOutput md::simulate_molecule(MoleculeOptions& options) {
    console::print_info("\nPreparing simulation for " + options.pdbfile.filename());

    //##################################//
    //###           GLOBALS          ###//
    //##################################//
    io::Folder mdp_folder = options.output + "mdp/";
    io::Folder setup_path = options.output + "protein/setup/";
    io::Folder em_path = options.output + "protein/em/";
    io::Folder eq_path = options.output + "protein/eq/";
    io::Folder prod_path = options.output + "protein/prod/";

    // gmx::gmx::set_outputlog(output + "gmx.log");
    //##################################//
    //###           SETUP            ###//
    //##################################//
    GROFile solv_ion(setup_path + "solv_ion.gro");
    TOPFile top(setup_path + "topol.top");
    if (!solv_ion.exists() || !top.exists()) {
        GROFile conf(setup_path + "conf.gro");
        if (!conf.exists()) {
            std::cout << "\tConverting structure to GROMACS format..." << std::flush;

            // prepare the pdb file for gromacs
            auto[conf, _, posre] = pdb2gmx(options.pdbfile)
                .output(setup_path)
                .ignore_hydrogens()
                // .virtual_sites()
                .water_model(options.watermodel)
                .forcefield(options.forcefield)
            .run();
            std::cout << " done." << std::endl;
        } else {
            std::cout << "\tReusing previously generated GROMACS structure." << std::endl;
        }

        // create a box around the protein
        std::cout << "\tGenerating unit cell..." << std::flush;
        auto[uc] = editconf(conf)
            .output(setup_path + "uc.gro")
            .box_type(options.boxtype)
            .extend(1.5)
        .run();
        std::cout << " done." << std::endl;

        // add water to the box
        std::cout << "\tSolvating unit cell..." << std::flush;
        auto[solv] = solvate(uc)
            .output(setup_path + "solv.gro")
            .solvent(options.forcefield, options.watermodel)
            .radius(0.2)
            .topology(top)
        .run();
        std::cout << " done." << std::endl;

        // generate an empty tpr file 
        MDPFile mdp = MDPFile(setup_path + "empty.mdp"); mdp.create();
        auto[ions] = grompp(mdp, top, solv)
            .output(setup_path + "ions.tpr")
            .warnings(1)
        .run();
        mdp.remove();

        // add ions to the box
        std::cout << "\tAdding ions to unit cell..." << std::flush;
        genion(ions)
            .output(solv_ion)
            .topology(top)
            .neutralize()
            .cation(options.cation)
            .anion(options.anion)
            .ion_group("SOL")
        .run();
        std::cout << " done." << std::endl;
    } else {
        std::cout << "\tReusing previously generated system setup." << std::endl;
    }

    //##################################//
    //###     ENERGY MINIMIZATION    ###//
    //##################################//
    GROFile emgro(em_path + "em.gro");
    if (!emgro.exists()) {
        std::cout << "\tRunning energy minimization..." << std::flush;

        // prepare energy minimization sim
        MDPFile mdp = EMMDPCreator().write(mdp_folder + "emmol.mdp");
        auto[emtpr] = grompp(mdp, top, solv_ion)
            .output(em_path + "em.tpr")
            .restraints(solv_ion)
        .run();

        // run energy minimization
        mdrun(emtpr)
            .output(em_path, "em")
            .jobname(options.name + "_mol")
        .run(options.setupsim, options.jobscript)->submit();
        std::cout << " done." << std::endl;
    } else {
        std::cout << "\tReusing previously generated energy minimization." << std::endl;
    }

    //##################################//
    //###       EQUILIBRATION        ###//
    //##################################//
    GROFile eqgro(eq_path + "eq.gro");
    NDXFile index(eq_path + "index.ndx");
    if (!eqgro.exists()) {
        std::cout << "\tRunning thermalization..." << std::flush;

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
        MDPFile mdp = EQMDPCreatorMol().write(mdp_folder + "eqmol.mdp");
        auto[eqtpr] = grompp(mdp, top, emgro)
            .output(eq_path + "eq.tpr")
            .index(index)
            .restraints(emgro)
            .warnings(1)
        .run();

        // run equilibration
        mdrun(eqtpr)
            .output(eq_path, "eq")
            .jobname(options.name + "_mol")
        .run(options.setupsim, options.jobscript)->submit();
        std::cout << " done." << std::endl;
    } else {
        std::cout << "\tReusing previously generated thermalization." << std::endl;
    }

    //##################################//
    //###       PRODUCTION           ###//
    //##################################//
    GROFile prodgro(prod_path + "confout.gro");
    if (!prodgro.exists()) {
        std::cout << "\tGenerating backbone restraints..." << std::flush;
        auto[backbone] = genrestr(eq_path + "eq.gro")
            .output(setup_path + "backbone.itp")
            .index(index)
            .force(2000, 2000, 2000)
        .run();
        auto backbone_chains = backbone.split_restraints(top.get_includes());
        top.include(backbone_chains, "POSRESBACKBONE");
        // top.include(top.relative_path(backbone), "POSRESBACKBONE", "chain topol");
        std::cout << " done." << std::endl;

        // prepare production sim
        std::cout << "\tRunning production..." << std::flush;
        MDPFile mdp = options.molmdp->write(mdp_folder + "prmol.mdp");
        auto[prodtpr] = grompp(mdp, top, eqgro)
            .output(prod_path + "prod.tpr")
            .index(index)
            .restraints(eqgro)
            .warnings(1)
        .run();

        // run production
        auto job = mdrun(prodtpr)
            .output(prod_path, "prod")
            .jobname(options.name + "_mol")
        .run(options.mainsim, options.jobscript);
        std::cout << " done." << std::endl;

        return {std::move(job), top, solv_ion};
    } else {
        std::cout << "\tReusing previous position-restrained simulation." << std::endl;
    }
    auto job = std::make_unique<NoExecution<MDRunResult>>(prod_path);
    return {std::move(job), top, solv_ion};
}