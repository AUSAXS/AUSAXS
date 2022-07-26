#include <CLI/CLI.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <data/Body.h>
#include <data/Protein.h>
#include <fitter/IntensityFitter.h>
#include <fitter/Fit.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>

using std::cout, std::endl;

int main(int argc, char const *argv[]) {
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};

    string input_structure, input_measurement, output = "figures/", placement_strategy;
    app.add_option("input_s", input_structure, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", input_measurement, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", setting::axes::scattering_intensity_plot_binned_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h", setting::grid::rh, "Radius of the hydration atoms.");
    app.add_option("--qlow", setting::fit::q_low, "Lower limit on used q values from measurement file.");
    app.add_option("--qhigh", setting::fit::q_high, "Upper limit on used q values from measurement file.");
    app.add_flag("--center,!--no-center", setting::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {setting::grid::psc = setting::grid::PlacementStrategyChoice::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::psc = setting::grid::PlacementStrategyChoice::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::psc = setting::grid::PlacementStrategyChoice::JanStrategy;}

    Protein protein(input_structure);
    protein.generate_new_hydration();
    hist::ScatteringHistogram h = protein.get_histogram();

    IntensityFitter fitter(input_measurement, h);
    std::shared_ptr<Fit> result = fitter.fit();

    // Fit plot
    plots::PlotIntensityFit plot_f(fitter);
    plot_f.save(output + "intensity_fit." + setting::figures::format);

    // Residual plot
    plots::PlotIntensityFitResiduals plot_r(fitter);
    plot_r.save(output + "residuals." + setting::figures::format);

    FitReporter::report(result);

    vector<double> q;
    for (double qv = 0; qv < 1; qv+=0.01) {q.push_back(qv);}

    SimpleDataset data = fitter.get_model_dataset(q);
    plots::PlotDataset plot_d(data);
    plot_d.save(output + "fitted_model." + setting::figures::format);

    // double intercept = result->params["b"];
    double I0 = fitter.get_intercept();
    double DrhoV2 = std::pow(protein.get_relative_charge(), 2);
    double re2 = pow(constants::radius::electron*constants::unit::cm, 2);
    double m = protein.get_absolute_mass()*constants::unit::mg;

    cout << "concentration is: " << I0*m/(DrhoV2*re2) << endl;
    return 0;
}