#include <plots/PlotIntensityFit.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>

plots::PlotIntensityFit::PlotIntensityFit(SimpleIntensityFitter& fitter) : Plot() {
    prepare_canvas();
    auto graphs = fitter.plot();
    plot(graphs);
}

plots::PlotIntensityFit::PlotIntensityFit(const Fit& fit) : Plot() {
    prepare_canvas();
    plot(fit.figures);
}

plots::PlotIntensityFit::PlotIntensityFit(const std::shared_ptr<Fit> fit) : Plot() {
    prepare_canvas();
    plot(fit->figures);
}

plots::PlotIntensityFit::~PlotIntensityFit() = default;

void plots::PlotIntensityFit::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotIntensityFit::plot(const Fit::Plots& graphs) const {
    PlotOptions options0, options1, options2;

    options2.set("markers", {{"color", kOrange+1}, {"markerstyle", 7}, {"title", "Fit"}, {"xlabel", "q"}, {"ylabel", "Intensity"}});
    options0.set("markers", {{"color", kBlack}, {"share_axis", true}});
    options1.set("line", {{"color", kBlack}, {"share_axis", true}});

    draw(graphs.data, options2, canvas);
    draw(graphs.intensity_interpolated, options0, canvas);
    draw(graphs.intensity, options1, canvas);
}

void plots::PlotIntensityFit::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
    canvas->SetLogy();
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
}