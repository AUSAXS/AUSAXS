// includes
#include <vector>
#include <string>
#include <boost/program_options.hpp>

#include "Protein.cpp"



int main(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description description("Usage: ./program <input path> <output path> <optionals>");
    description.add_options()
        ("help,h", "Show this message.")
        // ("input,i", po::value<string>()->required(), "Path to the input file.")
        // ("output,p", po::value<string>()->required(), "Path to the output file.")
        ("reduce,r", po::value<int>(), "The factor to reduce the number of generated water molecules by (default: 3).")
        ("width,w", po::value<double>(), "The distance between each grid point (default: 1). Lower widths increases the precision.");

    int reduce = 0;
    double width = 1;
    // string input, output;
    try {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
        po::notify(vm);

        // input = vm["input"].as<string>();
        // output = vm["output"].as<string>();
        if(vm.count("help")) {
            cout << description << endl;
            exit(0);
        }
        if (vm.count("reduce")) {
            reduce = vm["reduce"].as<int>();
            cout << "Reduction factor set to " << reduce << endl;
        }
        if (vm.count("width")) {
            width = vm["width"].as<double>();
            cout << "Width set to " << width << endl;
        }
    } catch ( const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // cout << input << endl;
    // cout << output << endl;
    // exit(1);

    Protein protein(argv[1]);
    protein.generate_new_hydration(reduce, width);
    protein.save(argv[2]);
    return 0;
}