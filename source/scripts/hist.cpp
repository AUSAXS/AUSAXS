// includes
#include <vector>
#include <string>
#include <boost/program_options.hpp>

#include "Protein.cpp"

#include "TH1D.h"
#include "TCanvas.h"

int reduce = 0;
double width = 1;
string input, output;

void parse_params(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description description("Usage: ./program <input path> <output path> <optionals>");
    description.add_options()
        ("help,h", "Show this message.")
        ("input,i", po::value<string>()->required(), "Path to the input file.")
        ("output,p", po::value<string>()->required(), "Path to the output file.")
        ("reduce,r", po::value<int>(), "The factor to reduce the number of generated water molecules by (default: 3).")
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
            reduce = vm["reduce"].as<int>();
            cout << "Reduction factor set to " << reduce << endl;
        }
        if (vm.count("width")) {
            width = vm["width"].as<double>();
            cout << "Width set to " << width << endl;
        }
    } catch ( const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char const *argv[]) {
    parse_params(argc, argv);
    Protein protein(argv[1]);
    protein.generate_new_hydration(reduce, width);
    Distances d = protein.calc_distances();

    unique_ptr<TCanvas> c = std::make_unique<TCanvas>("c", "canvas", 600, 600);
    unique_ptr<TH1D> h = std::make_unique<TH1D>("h", "hist", 100, -10, 10);
    for (int i = 0; i < d.pp.size(); i++) {
        h->Fill(d.pp[i], d.wpp[i]);
    }
    for (int i = 0; i < d.hh.size(); i++) {
        h->Fill(d.hh[i], d.whh[i]);
    }
    for (int i = 0; i < d.hp.size(); i++) {
        h->Fill(d.hp[i], d.whp[i]);
    }
    c->SetRightMargin(0.15);
    c->SaveAs(output.c_str());

    return 0;
}