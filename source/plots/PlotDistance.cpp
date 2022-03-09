#include "plots/PlotDistance.h"
#include "ScatteringHistogram.h"
#include "settings.h"

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>

using std::unique_ptr, std::string;

void plots::PlotDistance::save(string path) const {
    unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("PlotDistanceCanvas", "canvas", 600, 600);
    auto hists = d.plot_distance();
    
    // use some nicer colors
    hists[0]->SetLineColor(kOrange+1);
    hists[1]->SetLineColor(kAzure+1);
    hists[2]->SetLineColor(kGreen+1);
    hists[3]->SetLineColor(kBlack);

    // titles
    hists[3]->GetXaxis()->SetTitle("Distance [#AA]");
    hists[3]->GetXaxis()->CenterTitle();
    hists[3]->GetYaxis()->SetTitle("Count");
    hists[3]->GetYaxis()->CenterTitle();
    // hists[3]->GetYaxis()->SetTitleOffset();

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
    canvas->SetLeftMargin(0.15);
    canvas->SaveAs(path.c_str());
}