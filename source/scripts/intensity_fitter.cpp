// includes
#include <vector>
#include <string>
#include <iostream>

#include "data/Protein.h"
#include "fitter/IntensityFitter.cpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>

using std::cout, std::endl;

int main(int, char const *argv[]) {
    setting::grid::psc = setting::grid::RadialStrategy;
    // setting::grid::width = 0.5;
    // setting::grid::ra = 1.5;
    // setting::grid::rh = 1.5;

    Protein protein(argv[1]);
    protein.generate_new_hydration();
    std::shared_ptr<ScatteringHistogram> h = protein.get_distances();

    IntensityFitter fitter(argv[2], h);
    std::shared_ptr<Fitter::Fit> result = fitter.fit();

//*** FIT PLOT ***//
    std::unique_ptr<TCanvas> c1 = std::make_unique<TCanvas>("c1", "canvas", 600, 600);
    auto graphs = fitter.plot();

    // use some nicer colors
    graphs[0]->SetLineColor(kBlack);
    graphs[1]->SetMarkerColor(kBlack);
    graphs[2]->SetMarkerColor(kOrange+1);
    graphs[2]->SetLineColor(kOrange+1);

    graphs[0]->SetMarkerStyle(7);
    graphs[2]->SetMarkerStyle(7);

    // draw the graphs
    graphs[0]->Draw("AP"); // Axes points
    graphs[2]->Draw("SAME P"); // Point
    graphs[1]->Draw("SAME L"); // Line

    // setup the canvas and save the plot
    string path = string(argv[3]) + "intensity_fit.pdf";
    c1->SetTitle("Fit");
    c1->SetLogy();
    c1->SetLogx();
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    c1->SaveAs(path.c_str());

//*** RESIDUAL PLOT ***//
    std::unique_ptr<TCanvas> c2 = std::make_unique<TCanvas>("c2", "canvas", 600, 600);
    std::unique_ptr<TGraphErrors> graph = fitter.plot_residuals();
    std::unique_ptr<TLine> line = std::make_unique<TLine>(0, 0, 1, 0); // solid black line at x=0

    // use some nicer colors
    graph->SetMarkerColor(kOrange+1);
    graph->SetLineColor(kOrange+1);
    line->SetLineColor(kBlack);

    graph->SetMarkerStyle(7);

    // draw the graphs
    graph->Draw("AP");
    line->Draw("SAME");

    // setup the canvas and save the plot
    path = string(argv[3]) + "residuals.pdf";
    c2->SetTitle("Residuals");
    c2->SetLogx();
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());

    result->print();
    cout << "c is: " << result->params["a"]*protein.get_mass()/pow(constants::radius::electron, 2)*constants::unit::mg/pow(constants::unit::cm, 3) << endl;
    return 0;
}