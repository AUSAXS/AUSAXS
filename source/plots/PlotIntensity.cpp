#include <fitter/Fit.h>
#include <plots/PlotIntensity.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPad.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

plots::PlotIntensity::PlotIntensity(hist::ScatteringHistogram&& d, int color) : Plot(), d(std::move(d)) {
    prepare_canvas();
    initial_intensity_plot(color);
}

plots::PlotIntensity::PlotIntensity(const hist::ScatteringHistogram& d, int color) : Plot(), d(d) {
    prepare_canvas();
    initial_intensity_plot(color);
}

plots::PlotIntensity::PlotIntensity(const SAXSDataset& d) : Plot() {
    prepare_canvas();
    initial_intensity_plot(d);
}

plots::PlotIntensity::~PlotIntensity() = default;

void plots::PlotIntensity::initial_intensity_plot(int color) {
    if (d.q.empty()) {throw except::invalid_operation("Error in PlotIntensity::plot_intensity: Class was not initialized with a histogram, and it has not been manually set.");}

    // Debye scattering intensity
    std::shared_ptr<TH1D> hI_debye = d.plot_debye_scattering();
    PlotOptions options;
    options.xlabel = "q";
    options.ylabel = "Intensity";
    options.line_width = 3;
    options.color = color;

    ymin = hI_debye->GetMinimum();
    ymax = hI_debye->GetMaximum();
    hI_debye->SetAxisRange(ymin*0.9, ymax*1.1, "Y"); // fix the axis range so we can match it with the guinier approx
    
    draw(hI_debye, options, canvas);
}

void plots::PlotIntensity::initial_intensity_plot(const Dataset& data) {
    std::shared_ptr<TGraph> graph = data.plot();
    plots::PlotOptions options = data.plot_options;
    options.xlabel = "q";
    options.ylabel = "Intensity";

    ymin = graph->GetYaxis()->GetXmin();
    ymax = graph->GetYaxis()->GetXmax();
    graph->GetHistogram()->SetMinimum(ymin*0.9);
    graph->GetHistogram()->SetMaximum(ymax*1.1);
    draw(graph, options, canvas);
}

void plots::PlotIntensity::quick_plot(const hist::ScatteringHistogram& h, std::string path) {
    PlotIntensity plot(h);
    plot.save(path);
}

void plots::PlotIntensity::plot_intensity(const Dataset& data) {
    std::shared_ptr<TGraph> graph = data.plot();
    PlotOptions options(data.plot_options);
    options.use_existing_axes = true;

    linpad->cd();
    graph->GetHistogram()->SetMinimum(ymin*0.9);
    graph->GetHistogram()->SetMaximum(ymax*1.1);
    draw(graph, options, canvas);
}

void plots::PlotIntensity::plot_intensity(const std::shared_ptr<Fit> fit, const PlotOptions& plot_options) {
    std::shared_ptr<TGraph> graph = fit->figures.intensity.plot();
    PlotOptions options(plot_options);
    options.use_existing_axes = true;

    linpad->cd();
    graph->GetHistogram()->SetMinimum(ymin*0.9);
    graph->GetHistogram()->SetMaximum(ymax*1.1);
    draw(graph, options, canvas);
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
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotIntensity::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
    linpad = std::make_unique<TPad>("PlotIntensityPad1", "linpad", 0, 0, 1, 1); // create a drawing pad

    linpad->SetLeftMargin(0.19);
    linpad->SetLogx();
    linpad->SetLogy();

    linpad->Draw();
    linpad->cd();
}