#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/FitReporter.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Console.h>
#include <form_factor/DisplacedVolumeTable.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string s_pdb, s_mfile, s_settings;
    bool use_existing_hydration = false, fit_excluded_volume = false;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/fit_all_exv/")->group("General options");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads)->group("General options");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from the measurement file.")->default_val(settings::axes::qmax)->group("General options");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from the measurement file.")->default_val(settings::axes::qmin)->group("General options");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");

    // protein options group
    app.add_flag("--center,!--no-center", settings::molecule::center, "Decides whether the protein will be centered.")->default_val(settings::molecule::center)->group("Protein options");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::molecule::use_effective_charge, "Decides whether the effective atomic charge will be used.")->default_val(settings::molecule::use_effective_charge)->group("Protein options");
    app.add_flag("--use-existing-hydration,!--no-use-existing-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used.")->default_val(use_existing_hydration)->group("Protein options");
    app.add_flag("--fit-excluded-volume,!--no-fit-excluded-volume", fit_excluded_volume, "Decides whether the excluded volume will be fitted.")->default_val(fit_excluded_volume)->group("Protein options");

    // advanced options group
    app.add_option("--reduce,-r", settings::grid::water_scaling, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")->default_val(settings::grid::water_scaling)->group("Advanced options");
    app.add_option("--grid_width,--gw", settings::grid::width, "The distance between each grid point in Ångström. Lower widths increase the precision.")->default_val(settings::grid::width)->group("Advanced options");
    app.add_option_function<std::string>("--placement-strategy,--ps", [] (const std::string& s) {settings::detail::parse_option("placement_strategy", {s});}, "The placement strategy to use. Options: Radial, Axes, Jan.")->group("Advanced options");
    app.add_option("--exv_radius,--er", settings::grid::exv_radius, "The radius of the excluded volume sphere used for the grid-based excluded volume calculations in Ångström.")->default_val(settings::grid::exv_radius)->group("Advanced options");
    app.add_flag("--exit-on-unknown-atom,!--no-exit-on-unknown-atom", settings::molecule::throw_on_unknown_atom, "Decides whether the program will exit if an unknown atom is encountered.")->default_val(settings::molecule::throw_on_unknown_atom)->group("Advanced options");
    app.add_flag("--implicit-hydrogens,!--no-implicit-hydrogens", settings::molecule::implicit_hydrogens, "Decides whether implicit hydrogens will be added to the structure.")->default_val(settings::molecule::implicit_hydrogens)->group("Advanced options");
    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));

    //###################//
    //### PARSE INPUT ###//
    //###################//
    io::ExistingFile pdb(s_pdb), mfile(s_mfile), settings(s_settings);
    settings::general::output += mfile.stem() + "/";

    //######################//
    //### ACTUAL PROGRAM ###//
    //######################//
    std::vector<settings::hist::HistogramManagerChoice> loop;
    std::vector<std::string> loop_names;

    std::string volumes = "None";
    #if TRAUBE_FF
        loop = {
            settings::hist::HistogramManagerChoice::HistogramManagerMT,
            settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid,
            settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface,
            settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg, 
            settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit,
            settings::hist::HistogramManagerChoice::FoXSManager,
            settings::hist::HistogramManagerChoice::CrysolManager,
            settings::hist::HistogramManagerChoice::PepsiManager
        };

        loop_names = {
            "HistogramManagerMT",
            "HistogramManagerMTFFGrid",
            "HistogramManagerMTFFGridSurface",
            "HistogramManagerMTFFAvg", 
            "HistogramManagerMTFFExplicit",
            "FoXS",
            "CRYSOL",
            "Pepsi-SAXS"
        };
        volumes = "TRAUBE";
    #elif PONTIUS_FF
        loop = {settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit};
        loop_names = {"HistogramManagerMTFFExplicit", "FoXS"};
        volumes = "PONTIUS";
    #endif
    std::cout << volumes << std::endl;

    io::Folder(settings::general::output).create();
    std::ofstream out(settings::general::output + "/" + volumes + "_" + mfile.stem() + ".txt");
    bool printed_volume = false;
    auto perform_fit = [&] (const std::string& name, settings::hist::HistogramManagerChoice choice, bool fit_exv) {
        settings::hist::histogram_manager = choice;

        data::Molecule protein(pdb);
        protein.generate_new_hydration();
        if (!printed_volume) {out << "size: " << std::to_string(protein.size_atom()) << std::endl; printed_volume = true;}

        std::shared_ptr<fitter::HydrationFitter> fitter;
        if (fit_exv) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(mfile, protein.get_histogram());}
        else {fitter = std::make_shared<fitter::HydrationFitter>(mfile, protein.get_histogram());}
        auto result = fitter->fit();
        std::cout << name << ": " << result->fval/result->dof << std::endl;
        out << name << ": " << result->fval/result->dof << std::endl;
        fitter::FitReporter::save(result.get(), settings::general::output + volumes + "/" + name + ".txt");
    };

    for (unsigned int i = 0; i < loop.size(); ++i) {
        switch (loop[i]) {
            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                settings::molecule::use_effective_charge = true;
                perform_fit(loop_names[i], loop[i], false);
                break;
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid:
                settings::molecule::use_effective_charge = false;
                perform_fit(loop_names[i], loop[i], false);
                break;
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface:
                settings::molecule::use_effective_charge = false;
                perform_fit(loop_names[i], loop[i], false);
                settings::grid::rvol = 2.15;
                perform_fit(loop_names[i] + "_fitted_215", loop[i], true);
                settings::grid::rvol = 3.00;
                perform_fit(loop_names[i] + "_fitted_300", loop[i], true);
                break;
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
            case settings::hist::HistogramManagerChoice::FoXSManager:
            case settings::hist::HistogramManagerChoice::PepsiManager:
            case settings::hist::HistogramManagerChoice::CrysolManager:
                settings::molecule::use_effective_charge = false;
                perform_fit(loop_names[i], loop[i], false);
                perform_fit(loop_names[i] + "_fitted", loop[i], true);
                break;
            default:
                throw except::invalid_argument("Unknown histogram manager choice.");
        }
    }
    return 0;
}