#include <plots/PlotIntensityFitResiduals.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(SimpleIntensityFitter& fitter) : Plot() {
    prepare_canvas();
    SimpleDataset graph = fitter.plot_residuals();
    plot(graph);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const Fit& fit) : Plot() {
    prepare_canvas();
    plot(fit.residuals);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const std::shared_ptr<Fit> fit) : Plot() {
    prepare_canvas();
    plot(fit->residuals);
}

plots::PlotIntensityFitResiduals::~PlotIntensityFitResiduals() = default;

void plots::PlotIntensityFitResiduals::quick_plot(const std::shared_ptr<Fit> fit, std::string path) {
    PlotIntensityFitResiduals plot(fit);
    plot.save(path);
}

void plots::PlotIntensityFitResiduals::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotIntensityFitResiduals::plot(const SimpleDataset graph) const {
    std::unique_ptr<TLine> line = std::make_unique<TLine>(0, 0, graph.x().back(), 0); // solid black line at x=0
    PlotOptions options("markers", {{"color", kOrange+1}, {"markerstyle", 7}, {"title", "Residuals"}, {"xlabel", "q"}, {"ylabel", "Residual"}});

    draw(graph, options, canvas);
    line->SetLineColor(kBlack);
    line->DrawClone("SAME");
}

void plots::PlotIntensityFitResiduals::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
}