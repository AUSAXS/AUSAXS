// includes
#include <vector>
#include <string>
#include <iostream>
#include <boost/program_options.hpp>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>

#include "data/Body.h"
#include "data/Protein.h"
#include "fitter/IntensityFitter.h"
#include "plots/PlotIntensityFit.h"
#include "plots/PlotIntensityFitResiduals.h"

using std::cout, std::endl;

string input_structure, input_measurement, output;
void parse_params(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description description("Usage: ./program <input structure path> <input measurement path> <output path> <optionals>");
    description.add_options()
        ("help,h", "Show this message.")
        ("input_s", po::value<string>()->required(), "Path to the input structure file.")
        ("input_m", po::value<string>()->required(), "Path to the input measurement file.")
        ("output", po::value<string>()->required(), "Path to the output folder.")
        ("qlow", po::value<double>(), "The lower limit on the q-values used for the fit.")
        ("qhigh", po::value<double>(), "The upper limit on the q-values used for the fit.")
        ("center", po::value<bool>()->implicit_value(true), "Set true to center the structure, or false to not.")
        ("placement_strategy", po::value<string>(), "The placement strategy to use. Options: Radial, Axes, Jan.")
        ("width", po::value<double>(), "The grid bin width used internally in Ångström. Default: 1");

    // set positional parameters
    boost::program_options::positional_options_description pos_desc;
    pos_desc.add("input_s", 1); // input structure is always the first argument
    pos_desc.add("input_m", 1); // input measurement is always the second argument
    pos_desc.add("output", 2); // output is always the third argument
    po::command_line_parser parser{argc, argv};
    po::parsed_options parsed = parser.options(description).positional(pos_desc).run();
    
    try {
        po::variables_map vm;
        po::store(parsed, vm);
        po::notify(vm);

        input_structure = vm["input_s"].as<string>();
        input_measurement = vm["input_m"].as<string>();
        output = vm["output"].as<string>();
        if(vm.count("help")) {
            cout << description << endl;
            exit(0);
        }
        if (vm.count("qlow")) {
            setting::fit::q_low = vm["qlow"].as<double>();
            cout << "Lower limit on input q-values set to " << setting::fit::q_low << "." << endl;
        }
        if (vm.count("qhigh")) {
            setting::fit::q_high = vm["qhigh"].as<double>();
            cout << "Upper limit on input q-values set to " << setting::fit::q_high << "." << endl;
        }
        if (vm.count("center")) {
            setting::protein::center = vm["center"].as<bool>();
            cout << "Structure will " << (setting::protein::center ? "be centered." : "not be centered.") << endl;
        }
        if (vm.count("width")) {
            setting::grid::width = vm["width"].as<double>();
            cout << "Grid bin width set to " << setting::grid::width << "." << endl;
        }
        if (vm.count("placement_strategy")) {
            string parsed = vm["placement_strategy"].as<string>();
            if (parsed == "Radial") {
                cout << "Using radial placement strategy." << endl;
                setting::grid::psc = setting::grid::RadialStrategy;
            }
            else if (parsed == "Axes") {
                cout << "Using axes placement strategy." << endl;
                setting::grid::psc = setting::grid::AxesStrategy;
            }
            else if (parsed == "Jan") {
                cout << "Using Jan placement strategy." << endl;
                setting::grid::psc = setting::grid::JanStrategy;
            }
        }

    } catch (const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char const *argv[]) {
    parse_params(argc, argv);
    // setting::axes::scattering_intensity_plot_binned_width = 0.5;
    // setting::grid::ra = 1.5;
    // setting::grid::rh = 1.5;

    Protein protein(input_structure);
    protein.generate_new_hydration();
    std::shared_ptr<ScatteringHistogram> h = protein.get_histogram();

    IntensityFitter fitter(input_measurement, h);
    std::shared_ptr<Fitter::Fit> result = fitter.fit();

    // Fit plot
    PlotIntensityFit plot_f(fitter);
    plot_f.save(output + "intensity_fit." + setting::figures::format);

    // Residual plot
    PlotIntensityFitResiduals plot_r(fitter);
    plot_r.save(output + "residuals." + setting::figures::format);

    result->print();
    cout << "c is: " << result->params["a"]*protein.get_mass()/pow(constants::radius::electron, 2)*constants::unit::mg/pow(constants::unit::cm, 3) << endl;
    return 0;
}