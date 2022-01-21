#pragma once

#include "plots/Plot.h"
#include "ScatteringHistogram.h"
#include "settings.h"

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

class PlotIntensity : public Plot {
  public:
    PlotIntensity(const ScatteringHistogram& d) : Plot(), d(d) {}
    ~PlotIntensity() override = default;

    void save(const std::string& path) const override {
        unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
        
        // Debye scattering intensity
        canvas->cd();
        unique_ptr<TPad> linpad = std::make_unique<TPad>("linpad", "pad", 0, 0, 1, 1); // create a drawing pad
        linpad->Draw();
        linpad->SetLogx();
        linpad->SetLogy();
        linpad->cd();
        auto hI_debye = d.plot_debye_scattering();
        hI_debye->SetLineWidth(3);
        hI_debye->SetLineColor(kOrange+1);
        double ymin = hI_debye->GetMinimum();
        double ymax = hI_debye->GetMaximum();
        hI_debye->SetAxisRange(ymin*0.9, ymax*1.1, "Y"); // fix the axis range so we can match it with the guinier approx
        
        // titles
        linpad->SetLeftMargin(0.19);
        hI_debye->GetXaxis()->SetTitle("q");
        hI_debye->GetXaxis()->CenterTitle();
        hI_debye->GetYaxis()->SetTitle("Intensity");
        hI_debye->GetYaxis()->CenterTitle();
        hI_debye->GetYaxis()->SetTitleOffset(1.2);
        hI_debye->Draw("HIST L");

        // Guinier approximation
        // we have to create a second drawing pad since our scattering intensity is now log10 I(q)
        canvas->cd();
        unique_ptr<TPad> logpad = std::make_unique<TPad>("logpad", "pad", 0, 0, 1, 1); 
        logpad->Draw();
        logpad->SetFillStyle(4000); // make this second plot transparent (otherwise it'd overwrite the first one)
        logpad->SetFillColor(0);
        logpad->SetFrameFillStyle(4000);
        logpad->SetLogx();
        logpad->cd();
        logpad->SetLeftMargin(0.19);

        auto hI_guinier = d.plot_guinier_approx();
        double offset = log10(ymax) - hI_guinier->GetMaximum(); // the offset from the debye plot (free variable in the guinier approx)
        for (int i = 1; i < hI_guinier->GetNbinsX(); i++) {hI_guinier->SetBinContent(i, hI_guinier->GetBinContent(i)+offset);} // apply the offset
        hI_guinier->SetLineWidth(3);
        hI_guinier->SetLineColor(kAzure+1);
        hI_guinier->SetLineStyle(kDashed);
        hI_guinier->SetAxisRange(log10(ymin*0.9), log10(ymax*1.1), "Y"); // use same limits as the debye plot
        hI_guinier->SetNdivisions(3, "Y"); // use only 3 labels on the y axis
        hI_guinier->Draw("Y+ HIST L"); // Y+ creates a second axis on the right side

        // Vertical line at the Guinier gyration ratio
        double Rg = sqrt(d.calc_guinier_gyration_ratio_squared());
        unique_ptr<TLine> gyration_ratio = std::make_unique<TLine>(1./Rg, hI_guinier->GetMaximum(), 1./Rg, hI_guinier->GetMinimum());
        gyration_ratio->SetLineColor(kBlack);
        gyration_ratio->SetLineStyle(kDashed);
        gyration_ratio->SetLineWidth(3);
        gyration_ratio->Draw("SAME");

        // setup the canvas and save the plot
        canvas->SaveAs(path.c_str());
    }

  private:
    const ScatteringHistogram& d;
};