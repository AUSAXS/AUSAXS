// includes
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <string>

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Protein.h"
#include "plot_style.h"

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
            setting::protein::grid_width = vm["width"].as<double>();
            cout << "Width set to " << setting::protein::grid_width << endl;
        }
    } catch ( const std::exception& e ) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char const *argv[]) {
    parse_params(argc, argv);

    setting::grid::psc = setting::grid::RadialStrategy;

    Protein protein(argv[1]);
    protein.generate_new_hydration();
    Distances d = protein.calc_distances();

    setup_style();

// Distance plot
    const vector<int> axes = {600, 0, 60};
    unique_ptr<TCanvas> c1 = std::make_unique<TCanvas>("c1", "canvas", 600, 600);
    unique_ptr<TH1D> h_pp = std::make_unique<TH1D>("h_pp", "hist", axes[0], axes[1], axes[2]);
    unique_ptr<TH1D> h_hh = std::make_unique<TH1D>("h_hh", "hist", axes[0], axes[1], axes[2]);
    unique_ptr<TH1D> h_hp = std::make_unique<TH1D>("h_hp", "hist", axes[0], axes[1], axes[2]);
    unique_ptr<TH1D> h_tot = std::make_unique<TH1D>("h_tot", "hist", axes[0], axes[1], axes[2]);

    vector<double> p_pp(axes[0], 0);
    vector<double> p_hh(axes[0], 0);
    vector<double> p_hp(axes[0], 0);
    vector<double> p_tot(axes[0], 0);
    double width = (double) (axes[2]-axes[1])/axes[0]; // very important to cast this operation to a double - divison by two ints
    for (int i = 0; i < d.pp.size(); i++) {
        p_pp[std::round(d.pp[i]/width)] += d.wpp[i];
    }
    for (int i = 0; i < d.hh.size(); i++) {
        p_hh[std::round(d.hh[i]/width)] += d.whh[i];
    }
    for (int i = 0; i < d.hp.size(); i++) {
        p_hp[std::round(d.hp[i]/width)] += d.whp[i];
    }

    for (int i = 1; i < axes[0]; i++) {
        h_pp->SetBinContent(i, p_pp[i-1]);
        h_hh->SetBinContent(i, p_hh[i-1]);
        h_hp->SetBinContent(i, p_hp[i-1]);
        p_tot[i-1] = p_pp[i-1] + p_hh[i-1] + p_hp[i-1];
        h_tot->SetBinContent(i, p_tot[i-1]);
    }

    // use some nicer colors
    h_pp->SetLineColor(kOrange+1);
    h_hh->SetLineColor(kAzure+1);
    h_hp->SetLineColor(kGreen+1);
    h_tot->SetLineColor(kBlack);

    // draw the histograms on the canvas
    h_tot->Draw("HIST L");
    h_pp->Draw("SAME HIST L");
    h_hh->Draw("SAME HIST L");
    h_hp->Draw("SAME HIST L");

    // create a legend
    unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.6, 0.65, 0.9, 0.9);
    legend->AddEntry("h_tot", "Total", "l");
    legend->AddEntry("h_pp", "Atom-atom", "l");
    legend->AddEntry("h_hh", "Water-water", "l");
    legend->AddEntry("h_hp", "Atom-water", "l");
    legend->SetTextSize(0.04);
    legend->Draw();

    // setup the canvas and save the plot
    string path = output + "distances.pdf";
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    c1->SaveAs(path.c_str());

// Debye scattering intensity plot
    vector<double> Iq = protein.debye_scattering_intensity(axes, p_tot);
    unique_ptr<TCanvas> c2 = std::make_unique<TCanvas>("c2", "canvas", 600, 600);
    unique_ptr<TH1D> h_I = std::make_unique<TH1D>("h_I", "hist", Iq.size(), 0, 10);

    for (int i = 0; i < Iq.size(); i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h_I->SetBinContent(i+1, Iq[i]);
        // cout << "Bin " << i << ": " << Iq[i] << endl;
    }

    h_I->Draw("HIST L");

    // setup the canvas and save the plot
    path = output + "intensity.pdf";
    c2->SetLogx();
    c2->SetLogy();
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());

    return 0;
}