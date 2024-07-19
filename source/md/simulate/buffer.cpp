/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/simulate/buffer.h>
#include <md/programs/all.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/files/MDPCreator.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>

md::SimulateBufferOutput md::simulate_buffer(const BufferOptions& options) {
    //##################################//
    //###           GLOBALS          ###//
    //##################################//
    io::Folder mdp_folder = options.output + "mdp/";
    io::Folder setup_path = options.output + "buffer/setup/";
    io::Folder em_path = options.output + "buffer/em/";
    io::Folder eq_path = options.output + "buffer/eq/";
    io::Folder prod_path = options.output + "buffer/prod/";

    if (!options.refgro.exists()) {
        throw except::io_error("Unit cell reference file not found");
    }

    //##################################//
    //###           SETUP            ###//
    //##################################//
    GROFile gro(setup_path + "buffer.gro");
    TOPFile top(setup_path + "topol.top");
    if (!gro.exists() || !top.exists()) {
        {
            std::cout << "\tCopying unit cell from molecule..." << std::flush;
            std::string ff = option::to_string(options.forcefield) + ".ff";
            std::string wm = option::to_string(options.watermodel) + ".itp";
            top.create(
                "; Topology file for the buffer simulation.     \n"
                "; Include forcefield parameters                \n"
                "#include \"" + ff + "/forcefield.itp\"         \n"
                "                                               \n"
                "; Include water topology                       \n"
                "#include \"" + ff + "/" + wm + "\"             \n"
                "#ifdef POSRES_WATER                            \n"
                "; Position restraint for each water oxygen     \n"
                "[ position_restraints ]                        \n"
                ";  i funct       fcx        fcy        fcz     \n"
                "    1    1       1000       1000       1000    \n"
                "#endif                                         \n"
                "                                               \n"
                "; Include topology for ions                    \n"
                "#include \"" + ff + "/ions.itp\"               \n"
                "                                               \n"
                "[ system ]                                     \n"
                "; Name                                         \n"
                "BUFFER                                         \n"
                "                                               \n"
                "[ molecules ]                                  \n"
                "; Compound                                     \n"
            );
        }

        auto[ucgro] = editconf(options.refgro)
            .output(setup_path + "uc.gro")
            .box_type(options.boxtype)
            .extend(2)
        .run();

        // remove everything but the unit cell dimensions
        std::string unit_cell = ucgro.get_unit_cell();
        ucgro.remove();
        ucgro.create(
            " \n"
            "0\n" +
            unit_cell
        );
        std::cout << " done" << std::endl;

        // solvate the unit cell
        std::cout << "\tSolvating unit cell..." << std::flush;
        auto[solvatedgro] = solvate(ucgro)
            .output(gro)
            .solvent(options.forcefield, options.watermodel)
            .topology(top)
        .run();
        std::cout << " done" << std::endl;
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
        MDPFile mdp = EMMDPCreator().write(mdp_folder + "emsol.mdp");
        auto[emtpr] = grompp(mdp, top, gro)
            .output(em_path + "em.tpr")
            .restraints(gro)
            .warnings(1)
        .run();

        // run energy minimization
        mdrun(emtpr)
            .output(em_path, "em")
            .jobname(options.name + "_buf")
        .run(options.setupsim, options.jobscript)->wait();
        std::cout << " done" << std::endl;
    } else {
        std::cout << "\tReusing previously generated energy minimization." << std::endl;
    }

    //##################################//
    //###       EQUILIBRATION        ###//
    //##################################//
    GROFile eqgro(eq_path + "eq.gro");
    if (!eqgro.exists()) {
        eqgro = GROFile(eq_path + "confout.gro");
    }

    if (!eqgro.exists()) {
        std::cout << "\tRunning thermalization..." << std::flush;

        // prepare equilibration sim
        MDPFile mdp = EQMDPCreatorSol().write(mdp_folder + "eqsol.mdp");
        auto[eqtpr] = grompp(mdp, top, emgro)
            .output(eq_path + "eq.tpr")
            .restraints(emgro)
            .warnings(1)
        .run();

        // run equilibration
        mdrun(eqtpr)
            .output(eq_path, "eq")
            .jobname(options.name + "_buf")
        .run(options.mainsim, options.jobscript)->wait();
        std::cout << " done" << std::endl;
    } else {
        std::cout << "\tReusing previously generated thermalization." << std::endl;
    }

    //##################################//
    //###       PRODUCTION          ###//
    //##################################//
    GROFile prodgro(prod_path + "confout.gro");
    if (!prodgro.exists()) {
        std::cout << "\tRunning production..." << std::flush;

        // prepare production sim
        MDPFile mdp = options.bufmdp->write(mdp_folder + "prbuf.mdp");
        auto[prodtpr] = grompp(mdp, top, eqgro)
            .output(prod_path + "prod.tpr")
            .warnings(1)
        .run();

        // run production
        auto job = mdrun(prodtpr)
            .output(prod_path, "prod")
            .jobname(options.name + "_buf")
        .run(options.mainsim, options.jobscript);
        std::cout << " done." << std::endl;

        return {std::move(job), top, prodgro};
    } else {
        std::cout << "\tReusing previous free simulation." << std::endl;
    }
    auto job = std::make_unique<NoExecution<MDRunResult>>(prod_path);
    return {std::move(job), top, prodgro};
}