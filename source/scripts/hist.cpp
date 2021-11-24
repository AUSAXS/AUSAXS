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
    shared_ptr<Distances> d = protein.get_distances();

    setup_style();

// Distance plot
    unique_ptr<TCanvas> c1 = std::make_unique<TCanvas>("c1", "canvas", 600, 600);
    const vector<int> axes = {600, 0, 60};
    auto hists = d->plot_distance(axes);

    // use some nicer colors
    hists[0]->SetLineColor(kOrange+1);
    hists[1]->SetLineColor(kAzure+1);
    hists[2]->SetLineColor(kGreen+1);
    hists[3]->SetLineColor(kBlack);

    // draw the histograms on the canvas
    hists[3]->Draw("HIST L");
    hists[0]->Draw("SAME HIST L");
    hists[1]->Draw("SAME HIST L");
    hists[2]->Draw("SAME HIST L");

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
    unique_ptr<TCanvas> c2 = std::make_unique<TCanvas>("c2", "canvas", 600, 600);

    // Debye scattering intensity
    c2->cd();
    unique_ptr<TPad> linpad = std::make_unique<TPad>("linpad", "pad", 0, 0, 1, 1);
    linpad->Draw();
    linpad->SetLogx();
    linpad->SetLogy();
    linpad->cd();
    auto hI_debye = d->plot_debye_scattering();
    hI_debye->Draw("HIST L");

    // Guinier approximation
    c2->cd();
    unique_ptr<TPad> logpad = std::make_unique<TPad>("logpad", "pad", 0, 0, 1, 1);
    logpad->Draw();
    logpad->SetFillStyle(4000); // transparent
    logpad->SetFillColor(0);
    logpad->SetFrameFillStyle(4000);
    logpad->SetLogx();
    logpad->cd();
    auto hI_guinier = d->plot_guinier_approx();
    hI_guinier->Draw("Y+ HIST L");

    // Vertical line at the Guinier gyration ratio
    double Rg = sqrt(d->calc_guinier_gyration_ratio_squared());
    unique_ptr<TLine> gyration_ratio = std::make_unique<TLine>(1./Rg, hI_guinier->GetMaximum(), 1./Rg, hI_guinier->GetMinimum());
    gyration_ratio->SetLineColor(kRed);
    gyration_ratio->SetLineStyle(kDashed);
    gyration_ratio->SetLineWidth(3);
    gyration_ratio->Draw("SAME");

    // setup the canvas and save the plot
    path = output + "intensity.pdf";
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());

    return 0;
}