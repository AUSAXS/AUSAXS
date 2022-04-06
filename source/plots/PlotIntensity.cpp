#include <plots/PlotIntensity.h>
#include <Exceptions.h>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

plots::PlotIntensity::PlotIntensity(ScatteringHistogram&& d, int color) : Plot(), d(std::move(d)) {
    prepare_canvas();
    initial_intensity_plot(color);
}

plots::PlotIntensity::PlotIntensity(const ScatteringHistogram& d, int color) : Plot(), d(d) {
    prepare_canvas();
    initial_intensity_plot(color);
}

plots::PlotIntensity::PlotIntensity(const SAXSDataset& d, int color) : Plot() {
    prepare_canvas();
    initial_intensity_plot(d, color);
}

void plots::PlotIntensity::initial_intensity_plot(int color) {
    if (d.q.empty()) {throw except::invalid_operation("Error in PlotIntensity::plot_intensity: Class was not initialized with a histogram, and it has not been manually set.");}

    // Debye scattering intensity
    auto hI_debye = d.plot_debye_scattering();
    hI_debye->SetLineWidth(3);
    hI_debye->SetLineColor(color);
    ymin = hI_debye->GetMinimum();
    ymax = hI_debye->GetMaximum();
    hI_debye->SetAxisRange(ymin*0.9, ymax*1.1, "Y"); // fix the axis range so we can match it with the guinier approx
    
    // titles
    hI_debye->GetXaxis()->SetTitle("q");
    hI_debye->GetXaxis()->CenterTitle();
    hI_debye->GetYaxis()->SetTitle("Intensity");
    hI_debye->GetYaxis()->CenterTitle();
    hI_debye->GetYaxis()->SetTitleOffset(1.2);
    hI_debye->DrawClone("HIST L");
}

void plots::PlotIntensity::initial_intensity_plot(const Dataset& data, int color) {
    auto graph = data.plot();
    graph->SetLineWidth(3);
    graph->SetLineColor(color);
    ymin = graph->GetYaxis()->GetXmin();
    ymax = graph->GetYaxis()->GetXmax();
    graph->GetYaxis()->SetRangeUser(ymin*0.9, ymax*1.1);
    
    // titles
    graph->GetXaxis()->SetTitle("q");
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->SetTitle("Intensity");
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->DrawClone("AP");
}

void plots::PlotIntensity::plot_intensity(const Dataset& data, int color, double alpha) {
    auto graph = data.plot();

    if (data.draw_as_line) {
        graph->SetLineColorAlpha(color, alpha);
        graph->DrawClone("SAME L");
    } else {
        graph->SetMarkerStyle(7);
        graph->SetMarkerColorAlpha(color, alpha);
        graph->DrawClone("SAME P");
    }
}

void plots::PlotIntensity::plot_guinier_approx() {
    if (d.q.empty()) {throw except::invalid_operation("Error in PlotIntensity::plot_guinier_approx: Class was not initialized with a histogram, and it has not been manually set.");}

    // Guinier approximation
    // we have to create a second drawing pad since our scattering intensity is now log10 I(q)
    canvas->cd();
    logpad = std::make_unique<TPad>("PlotIntensityPad2", "logpad", 0, 0, 1, 1); 
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
    hI_guinier->DrawClone("Y+ HIST L"); // Y+ creates a second axis on the right side

    // Vertical line at the Guinier gyration ratio
    double Rg = sqrt(d.calc_guinier_gyration_ratio_squared());
    unique_ptr<TLine> gyration_ratio = std::make_unique<TLine>(1./Rg, hI_guinier->GetMaximum(), 1./Rg, hI_guinier->GetMinimum());
    gyration_ratio->SetLineColor(kBlack);
    gyration_ratio->SetLineStyle(kDashed);
    gyration_ratio->SetLineWidth(3);
    gyration_ratio->DrawClone("SAME");
}

void plots::PlotIntensity::save(std::string path) const {
    canvas->SaveAs(path.c_str());
}

void plots::PlotIntensity::prepare_canvas() {
    canvas = std::make_unique<TCanvas>("PlotIntensityCanvas", "canvas", 600, 600);
    linpad = std::make_unique<TPad>("PlotIntensityPad1", "linpad", 0, 0, 1, 1); // create a drawing pad

    linpad->SetLeftMargin(0.19);
    linpad->SetLogx();
    linpad->SetLogy();

    linpad->Draw();
    linpad->cd();
}