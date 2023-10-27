#include <CLI/CLI.hpp>

#include <em/ImageStack.h>
#include <data/Body.h>
#include <data/Molecule.h>
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
#include <data/record/Atom.h>
#include <data/Body.h>
#include <data/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>
#include <settings/MoleculeSettings.h>
#include <fitter/HydrationFitter.h>
#include <form_factor/FormFactor.h>
#include <hydrate/Grid.h>

#include <cassert>

//************************************************************************************************* 
//********************************* PLOT SCATTERING STUFF *****************************************
//*************************************************************************************************
int main(int argc, char const *argv[]) {
    settings::axes::qmin = 0.02;
    settings::axes::qmax = 1;

    std::string base_path = "temp/debug/";
    SimpleDataset I_aa(base_path + "ausaxs/ausaxs_aa.dat");
    SimpleDataset I_ax(base_path + "ausaxs/ausaxs_ax.dat");
    SimpleDataset I_xx(base_path + "ausaxs/ausaxs_xx.dat");
    SimpleDataset I_aw(base_path + "ausaxs/ausaxs_aw.dat");
    SimpleDataset I_wx(base_path + "ausaxs/ausaxs_wx.dat");
    SimpleDataset I_ww(base_path + "ausaxs/ausaxs_ww.dat");

    SimpleDataset coords(base_path + "COORDINATE.dat");
    SimpleDataset C_xx(base_path + "exclvol.dat");
    SimpleDataset C_aa(base_path + "vaccum.dat");
    SimpleDataset foxs(base_path + "foxs_vacuum.dat");

    SimpleDataset foxs_aa(base_path + "foxs_aa.dat");
    SimpleDataset foxs_ax(base_path + "foxs_ax.dat");
    SimpleDataset foxs_xx(base_path + "foxs_xx.dat");
    SimpleDataset foxs_aw(base_path + "foxs_aw.dat");
    SimpleDataset foxs_wx(base_path + "foxs_wx.dat");
    SimpleDataset foxs_ww(base_path + "foxs_ww.dat");
    
    // SimpleDataset lyshr(base_path + "LYSHR.RSR");

    I_aa.normalize(1);
    I_ax.normalize(1);
    I_xx.normalize(1);
    coords.normalize(1);
    C_xx.normalize(1);
    C_aa.normalize(1);
    foxs.normalize(1);

    SimpleDataset I_sum(I_aa);
    for (unsigned int i = 0; i < I_aa.size(); ++i) {
        I_aa.y(i) = std::abs(I_aa.y(i));
        I_ax.y(i) = std::abs(I_ax.y(i));
        I_xx.y(i) = std::abs(I_xx.y(i));
        I_aw.y(i) = std::abs(I_aw.y(i));
        I_wx.y(i) = std::abs(I_wx.y(i));
        I_ww.y(i) = std::abs(I_ww.y(i));

        double cy = coords.interpolate_y(I_aa.x(i));
        I_sum.y(i) = I_aa.y(i) + I_ax.y(i) + I_xx.y(i);
        // I_aa.y(i) /= cy;
        // I_ax.y(i) /= cy;
        // I_xx.y(i) /= cy;
    }

    for (unsigned int i = 0; i < foxs_aa.size(); ++i) {
        foxs_aa.y(i) = std::abs(foxs_aa.y(i));
        foxs_ax.y(i) = std::abs(foxs_ax.y(i));
        foxs_xx.y(i) = std::abs(foxs_xx.y(i));
        foxs_aw.y(i) = std::abs(foxs_aw.y(i));
        foxs_wx.y(i) = std::abs(foxs_wx.y(i));
        foxs_ww.y(i) = std::abs(foxs_ww.y(i));
    }

    I_aa.normalize(1);
    I_ax.normalize(1);
    I_xx.normalize(1);
    I_aw.normalize(1);
    I_wx.normalize(1);
    I_ww.normalize(1);

    foxs_aa.normalize(1);
    foxs_ax.normalize(1);
    foxs_xx.normalize(1);
    foxs_aw.normalize(1);
    foxs_wx.normalize(1);
    foxs_ww.normalize(1);

    coords.normalize(1);
    C_xx.normalize(1);
    C_aa.normalize(1);
    I_sum.normalize(1);
    foxs.normalize(1);

    I_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}});
    I_ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ax"}, {"color", style::color::blue}});
    I_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::green}});
    I_aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aw"}, {"color", style::color::pink}});
    I_wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xw"}, {"color", style::color::purple}});
    I_ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ww"}, {"color", style::color::brown}});

    foxs_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aa"}, {"color", style::color::red}});
    foxs_ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ax"}, {"color", style::color::blue}});
    foxs_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xx"}, {"color", style::color::green}});
    foxs_aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aw"}, {"color", style::color::pink}});
    foxs_wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xw"}, {"color", style::color::purple}});
    foxs_ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ww"}, {"color", style::color::brown}});

    coords.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "coords"}, {"linestyle", style::line::dashed}, {"color", style::color::red}});
    C_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_xx"}, {"linestyle", style::line::dashed}, {"color", style::color::blue}});
    C_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_aa"}, {"linestyle", style::line::dashed}, {"color", style::color::green}});
    I_sum.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "I_sum"}, {"color", style::color::black}});
    foxs.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "foxs"}, {"color", style::color::black}});
    // lyshr.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "lyshr"}, {"color", style::color::black}});

    // I_aa.save(base_path + "AUSAXS_aa.dat");
    // I_ax.save(base_path + "AUSAXS_ax.dat");
    // I_xx.save(base_path + "AUSAXS_xx.dat");
    // exclvol.save(base_path + "scaled_exclvol.dat");
    // vacuum.save(base_path + "scaled_vacuum.dat");

    plots::PlotDataset(I_aa)
        .plot(I_xx)
        .plot(C_xx)
        .plot(C_aa)
    .save(base_path + "compare.png");

    plots::PlotDataset(I_aa)
        .plot(C_aa)
        .plot(foxs)
    .save(base_path + "vacuum.png");

    plots::PlotDataset(I_aa)
        .plot(I_xx)
        .plot(I_ww)
        .plot(foxs_aa)
        .plot(foxs_xx)
        .plot(foxs_ww)
    .save(base_path + "compare_foxs.png");

    plots::PlotDataset(I_ax)
        .plot(I_aw)
        .plot(I_wx)
        .plot(foxs_ax)
        .plot(foxs_aw)
        .plot(foxs_wx)
    .save(base_path + "compare_foxs_cross.png");

    plots::PlotDataset(foxs_xx)
        .plot(C_xx)
        .plot(I_xx)
    .save(base_path + "exv.png");

    for (unsigned int i = 0; i < C_xx.size(); ++i) {
        C_aa.y(i) /= I_aa.interpolate_y(C_aa.x(i));
        C_xx.y(i) /= I_xx.interpolate_y(C_xx.x(i));
    }
    C_xx.normalize();
    C_aa.normalize();

    for (unsigned int i = 0; i < foxs_xx.size(); ++i) {
        foxs_aa.y(i) /= I_aa.interpolate_y(foxs_aa.x(i));
        foxs_xx.y(i) /= I_xx.interpolate_y(foxs_xx.x(i));
    }
    foxs_aa.normalize();
    foxs_xx.normalize();

    foxs_xx.add_plot_options({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}});
    plots::PlotDataset(foxs_xx)
        .plot(C_xx)
        .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(base_path + "exv_normalized.png");

    foxs_aa.add_plot_options({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}});
    plots::PlotDataset(foxs_aa)
        .plot(C_aa)
        .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(base_path + "vacuum_normalized.png");


    // for (unsigned int i = 0; i < foxs_xx.size(); ++i) {foxs_xx.y(i) /= foxs_aa.y(i);}
    // for (unsigned int i = 0; i < I_xx.size(); ++i) {I_xx.y(i) /= I_aa.y(i);}

    // plots::PlotDataset(foxs_xx)
    //     .plot(I_xx)
    //     .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    // .save(base_path + "xx_div_aa.png");
}

//*************************************************************************************************
//********************************* PLOT FORM FACTORS *********************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::general::output = "temp/stuff/ff/";
//     auto plot = [] (form_factor::FormFactor ff, std::string name) {
//         auto q = Axis(0, 1, 100).as_vector();
//         std::vector<double> y(q.size());
//         for (unsigned int i = 0; i < q.size(); ++i) {
//             y[i] = ff.evaluate(q[i]);
//         }
//         SimpleDataset d(q, y);
//         d.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff"}});
//         return d;
//     };

//     auto C = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::C), "carbon");
//     auto O = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::O), "oxygen");
//     auto EV = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::EXCLUDED_VOLUME), "excluded_volume");
//     auto other = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::OTHER), "other");
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

//*************************************************************************************************
//********************************* PLOT EXCLUDED VOLUME ******************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::protein::use_effective_charge = false;
//     io::ExistingFile file(argv[1]);
//     Protein protein(file);
//     protein.get_grid()->expand_volume();
//     protein.get_grid()->save("temp/stuff/grid.pdb");
// }

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

//*************************************************************************************************
//********************************* DEBUG SEQUENCER ***********************************************
//*************************************************************************************************
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