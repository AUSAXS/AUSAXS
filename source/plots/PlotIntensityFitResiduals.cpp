#include <plots/PlotIntensityFitResiduals.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(SimpleIntensityFitter& fitter) : Plot() {
    prepare_canvas();
    std::shared_ptr<TGraph> graph = fitter.plot_residuals();
    plot(graph);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const Fit& fit) : Plot() {
    prepare_canvas();
    plot(fit.residual_plot);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const std::shared_ptr<Fit> fit) : Plot() {
    prepare_canvas();
    plot(fit->residual_plot);
}

void plots::PlotIntensityFitResiduals::save(std::string path) const {
    canvas->SaveAs(path.c_str());
}

void plots::PlotIntensityFitResiduals::plot(const std::shared_ptr<TGraph> graph) const {
    std::unique_ptr<TLine> line = std::make_unique<TLine>(0, 0, graph->GetXaxis()->GetXmax(), 0); // solid black line at x=0

    // use some nicer colors
    graph->SetMarkerColor(kOrange+1);
    graph->SetLineColor(kOrange+1);
    line->SetLineColor(kBlack);

    graph->SetMarkerStyle(7);

    // set titles
    graph->SetTitle("Residuals");
    graph->GetXaxis()->SetTitle("q");
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleOffset(1.05);
    graph->GetYaxis()->SetTitle("Residual");
    graph->GetYaxis()->CenterTitle();

    // draw the graphs
    graph->DrawClone("AP");
    line->DrawClone("SAME");
}

void plots::PlotIntensityFitResiduals::prepare_canvas() {
    canvas = std::make_unique<TCanvas>("PlotIntensityFitResidualsCanvas", "canvas", 600, 600);
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
}