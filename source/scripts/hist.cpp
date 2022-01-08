// includes
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <string>

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"

#include "data/Body.h"
#include "data/Protein.h"
#include "plot_style.h"
#include "plots/PlotDistance.h"
#include "plots/PlotIntensity.h"

string input, output;
void parse_params(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description description("Usage: ./program <input path> <output path> <optionals>");
    description.add_options()
        ("help,h", "Show this message.")
        ("input,i", po::value<string>()->required(), "Path to the input file.")
        ("output,p", po::value<string>()->required(), "Path to the output file.")
        ("placement_strategy", po::value<string>(), "The placement strategy to use. Options: Radial, Axes, Jan.")
        ("reduce,r", po::value<double>(), "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")
        ("width,w", po::value<double>(), "The distance between each grid point (default: 1). Lower widths increases the precision.");

    // set positional parameters
    boost::program_options::positional_options_description pos_desc;
    pos_desc.add("input", 1); // input is always the first argument
    pos_desc.add("output", 2); // output is always the second argument
    po::command_line_parser parser{argc, argv};
    po::parsed_options parsed = parser.options(description).positional(pos_desc).run();
    
    try {
        po::variables_map vm; 
        po::store(parsed, vm);
        po::notify(vm);

        input = vm["input"].as<string>();
        output = vm["output"].as<string>();
        if(vm.count("help")) {
            cout << description << endl;
            exit(0);
        }
        if (vm.count("reduce")) {
            setting::grid::percent_water = vm["reduce"].as<int>();
            cout << "Percentage of water set to " << setting::grid::percent_water << endl;
        }
        if (vm.count("width")) {
            setting::grid::width = vm["width"].as<double>();
            cout << "Width set to " << setting::grid::width << endl;
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

    } catch ( const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char const *argv[]) { 
    parse_params(argc, argv);

    // setting::axes::scattering_intensity_plot_binned_width = 0.5;
    // setting::figures::format = "png";

    Protein protein(input);
    protein.generate_new_hydration();
    shared_ptr<ScatteringHistogram> d = protein.get_histogram();

    // Distance plot
    PlotDistance d_plot(d);
    d_plot.save(output + "distances." + setting::figures::format); 

    // Debye scattering intensity plot
    PlotIntensity i_plot(d);
    i_plot.save(output + "intensity." + setting::figures::format);
    return 0;
}