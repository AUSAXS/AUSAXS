#include <CLI/CLI.hpp>

#include <em/ImageStack.h>
#include <data/Body.h>
#include <data/Protein.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/File.h>
#include <rigidbody/sequencer/Sequencer.h>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <data/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>
#include <settings/ProteinSettings.h>
#include <fitter/HydrationFitter.h>
#include <hist/detail/FormFactor.h>
#include <hydrate/Grid.h>

#include <cassert>

// int main(int argc, char const *argv[]) {
//     settings::general::output = "temp/stuff/ff/";
//     auto plot = [] (hist::detail::FormFactor ff, std::string name) {
//         auto q = Axis(0, 1, 100).as_vector();
//         std::vector<double> y(q.size());
//         for (unsigned int i = 0; i < q.size(); ++i) {
//             y[i] = ff.evaluate(q[i]);
//         }
//         SimpleDataset d(q, y);
//         d.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff"}});
//         return d;
//     };

//     auto C = plot(hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::NEUTRAL_CARBON), "carbon");
//     auto O = plot(hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::NEUTRAL_OXYGEN), "oxygen");
//     auto EV = plot(hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::EXCLUDED_VOLUME), "excluded_volume");
//     auto other = plot(hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::OTHER), "other");
//     C.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff"}, {"legend", "Carbon"}, {"color", "red"}});
//     O.add_plot_options({{"legend", "Oxygen"}, {"color", "blue"}});
//     EV.add_plot_options({{"legend", "Excluded volume"}, {"color", "green"}});
//     other.add_plot_options({{"legend", "Other"}, {"color", "black"}});

//     plots::PlotDataset(C)
//         .plot(O)
//         .plot(EV)
//         .plot(other)
//     .save(settings::general::output + "ff.png");
// }

// int main(int argc, char const *argv[]) {
//     settings::protein::use_effective_charge = false;
//     io::ExistingFile mapfile(argv[1]);
//     io::ExistingFile mfile(argv[2]);
//     em::ImageStack images(mapfile);

//     settings::em::alpha_levels = {14, 15};
//     settings::fit::max_iterations = 100;
//     Axis axis(settings::em::alpha_levels, settings::fit::max_iterations);
//     fitter::HydrationFitter fitter(mfile);
//     unsigned int c = 0;
//     for (auto& level : axis.as_vector()) {
//         auto protein = images.get_protein(level);
//         protein->generate_new_hydration();
//         fitter.set_scattering_hist(protein->get_histogram());
//         auto res = fitter.fit();
//         std::cout << "Step " << ++c << ": Evaluated cutoff value " << level << " with chi2 " << res->fval << std::endl;
//     }
// }

int main(int argc, char const *argv[]) {
    settings::protein::use_effective_charge = false;
    io::ExistingFile file(argv[1]);
    Protein protein(file);
    protein.get_grid()->expand_volume();
    protein.get_grid()->save("temp/stuff/grid.pdb");
}

// int main(int argc, char const *argv[]) {
//     settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
//     settings::protein::use_effective_charge = false;

//         int distance = settings::rigidbody::bond_distance;
//         Atom a1 = Atom(Vector3<double>(0, 0, 0*distance), 1, "C", "C", 1);
//         Atom a2 = Atom(Vector3<double>(0, 0, 1*distance), 1, "C", "C", 1);
//         Atom a3 = Atom(Vector3<double>(0, 0, 2*distance), 1, "C", "C", 1);
//         Atom a4 = Atom(Vector3<double>(0, 0, 3*distance), 1, "C", "C", 1);

//         Body b1 = Body(std::vector<Atom>{a1});
//         Body b2 = Body(std::vector<Atom>{a2});
//         Body b3 = Body(std::vector<Atom>{a3});
//         Body b4 = Body(std::vector<Atom>{a4});
//         std::vector<Body> ap = {b1, b2, b3, b4};
//         rigidbody::RigidBody rigidbody(ap);
//         assert(rigidbody.get_constraint_manager()->distance_constraints.size() == 3);

//         rigidbody::RigidBody rigidbody2 = rigidbody::BodySplitter::split("data/LAR1-2/LAR1-2.pdb", {9, 99});
//         assert(rigidbody2.get_constraint_manager()->distance_constraints.size() == 2);
// }

// int main(int argc, char const *argv[]) {
//     settings::grid::scaling = 2;
//     settings::grid::cubic = true;
//     settings::general::verbose = true;
//     settings::axes::distance_bin_width = 0.1;

//     CLI::App app{"Rigid-body optimization."};
//     io::File pdb, mfile, settings;
//     std::vector<unsigned int> constraints;
//     app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
//     app.add_option("input_m", mfile, "Path to the measuremed data.")->required()->check(CLI::ExistingFile);
//     app.add_option("output", settings::general::output, "Path to save the hydrated file at.")->default_val("output/rigidbody/");
//     auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile);
//     app.add_option("--reduce,-r", settings::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
//     app.add_option("--grid_width,-w", settings::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
//     app.add_option("--bin_width", settings::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
//     app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
//     app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
//     auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
//     app.add_option("--iterations", settings::rigidbody::iterations, "Maximum number of iterations. Default: 1000.");
//     app.add_option("--constraints", settings::rigidbody::detail::constraints, "Constraints to apply to the rigid body.");
//     app.add_flag("--center,!--no-center", settings::protein::center, "Decides whether the protein will be centered. Default: true.");
//     app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
//     CLI11_PARSE(app, argc, argv);
    
//     //###################//
//     //### PARSE INPUT ###//
//     //###################//
//     settings::general::output += mfile.stem() + "/";

//     // if a settings file was provided
//     if (p_settings->count() != 0) {
//         settings::read(settings);        // read it
//         CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
//     } else {                            // otherwise check if there is a settings file in the same directory
//         if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
//             CLI11_PARSE(app, argc, argv);
//         }
//     }
//     if (settings::rigidbody::detail::constraints.empty()) {
//         throw except::missing_option("rigidbody: Constraints must be supplied. Use --constraints to specify them.");
//     }

//     rigidbody::sequencer::Sequencer(mfile, rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints))
//         .body_select_strategy(settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect)
//         .parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly)
//             .decay_strategy(settings::rigidbody::DecayStrategyChoice::Exponential)
//             .amplitude(0.5)
//         .transform_strategy(settings::rigidbody::TransformationStrategyChoice::RigidTransform)
//         .loop(5)
//     .execute();

//     return 0;
// }