#pragma once

#include "plots/Plot.h"
#include "ScatteringHistogram.h"
#include "fitter/IntensityFitter.h"
#include "settings.h"

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

class PlotIntensityFitResiduals : public Plot {
  public:
    PlotIntensityFitResiduals(IntensityFitter& fitter) : Plot(), fitter(fitter) {}
    ~PlotIntensityFitResiduals() override = default;

    void save(const std::string& path) const override {
        std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
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

  private:
    IntensityFitter& fitter;
};