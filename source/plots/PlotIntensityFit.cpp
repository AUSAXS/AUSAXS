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

void plots::PlotIntensityFit::plot(const Multiset& graphs) const {
    PlotOptions options0, options1, options2;

    options2.set("markers", {{"color", kOrange+1}, {"markerstyle", 7}, {"title", "Fit"}, {"xlabel", "q"}, {"ylabel", "Intensity"}});
    options1.set("line", {{"color", kBlack}, {"share_axis", true}});
    options0.set("markers", {{"color", kBlack}, {"share_axis", true}});

    draw(graphs[2], options2);
    draw(graphs[0], options0);
    draw(graphs[1], options1);
}

void plots::PlotIntensityFit::prepare_canvas() {
    canvas = std::make_unique<TCanvas>("PlotIntensityFitCanvas", "canvas", 600, 600);
    canvas->SetLogy();
    canvas->SetLogx();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);
}