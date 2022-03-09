#pragma once

#include <plots/PlotIntensityFitResiduals.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

void plots::PlotIntensityFitResiduals::save(std::string path) const {
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("PlotIntensityFitResidualsCanvas", "canvas", 600, 600);
    std::unique_ptr<TGraphErrors> graph = fitter.plot_residuals();
    std::unique_ptr<TLine> line = std::make_unique<TLine>(0, 0, graph->GetXaxis()->GetXmax(), 0); // solid black line at x=0

    // use some nicer colors
    graph->SetMarkerColor(kOrange+1);
    graph->SetLineColor(kOrange+1);
    line->SetLineColor(kBlack);

    graph->SetMarkerStyle(7);

    // draw the graphs
    graph->Draw("AP");
    line->Draw("SAME");

    // set titles
    graph->SetTitle("Residuals");
    graph->GetXaxis()->SetTitle("q");
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleOffset(1.05);
    graph->GetYaxis()->SetTitle("Residual");
    graph->GetYaxis()->CenterTitle();

    // setup the canvas and save the plot
    canvas->SetTitle("Residuals");
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
    canvas->SaveAs(path.c_str());
}