#include "plots/PlotIntensityFit.h"

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

void plots::PlotIntensityFit::save(std::string path) const {
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("PlotIntensityFitCanvas", "canvas", 600, 600);
    auto graphs = fitter.plot();

    // use some nicer colors
    graphs[0]->SetLineColor(kBlack);
    graphs[1]->SetMarkerColor(kBlack);
    graphs[2]->SetMarkerColor(kOrange+1);
    graphs[2]->SetLineColor(kOrange+1);

    graphs[0]->SetMarkerStyle(7);
    graphs[2]->SetMarkerStyle(7);

    // set titles
    graphs[0]->SetTitle("Fit");
    graphs[0]->GetXaxis()->SetTitle("q");
    graphs[0]->GetXaxis()->CenterTitle();
    graphs[0]->GetXaxis()->SetTitleOffset(1.05);
    graphs[0]->GetYaxis()->SetTitle("Intensity");
    graphs[0]->GetYaxis()->CenterTitle();

    // draw the graphs
    graphs[2]->Draw("AP"); // Point
    graphs[0]->Draw("SAME P"); // Axes points
    graphs[1]->Draw("SAME L"); // Line

    // setup the canvas and save the plot
    canvas->SetLogy();
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
    canvas->SaveAs(path.c_str());
}