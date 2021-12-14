// includes
#include <vector>
#include <string>
#include <boost/program_options.hpp>

#include "data/Protein.h"
#include "settings.h"

string input, output;
void parse_params(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description description("Usage: ./program <input path> <output path> <optionals>");
    description.add_options()
        ("help,h", "Show this message.")
        ("input,i", po::value<string>()->required(), "Path to the input file.")
        ("output,p", po::value<string>()->required(), "Path to the output file.")
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
    } catch (const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char const *argv[]) {
    parse_params(argc, argv);

    setting::grid::rh = 3;
    // setting::grid::psc = setting::grid::RadialStrategy;
    // setting::grid::percent_water = 0;

    Protein protein(argv[1]);
    protein.generate_new_hydration();
    protein.save(argv[2]);
    return 0;
}