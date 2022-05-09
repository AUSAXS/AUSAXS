#include <plots/PlotHistogram.h>
#include <utility/Utility.h>

#include <TCanvas.h>
#include <TH1D.h>

plots::PlotHistogram::PlotHistogram(const hist::Histogram& h) {
    prepare_canvas();
    initial_plot(h);
}

plots::PlotHistogram::~PlotHistogram() = default;

void plots::PlotHistogram::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotHistogram::quick_plot(const hist::Histogram& hist, std::string path) {
    plots::PlotHistogram plot(hist);
    plot.save(path);
}

void plots::PlotHistogram::initial_plot(const hist::Histogram& hist) {
    plots::PlotOptions options;
    draw(hist, options, canvas);
}

void plots::PlotHistogram::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
}