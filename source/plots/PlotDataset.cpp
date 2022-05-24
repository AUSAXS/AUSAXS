#include <plots/PlotDataset.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TCanvas.h>

using std::string, std::vector;

plots::PlotDataset::PlotDataset(const Dataset& data) : Plot() {
    prepare_canvas();
    initial_plot(data);
}

plots::PlotDataset::PlotDataset(const Multiset& data) {
    if (data.empty()) {throw except::size_error("Error in PlotDataset::PlotDataset: The given Multiset is empty!");}

    prepare_canvas();
    draw(data, canvas);
}

plots::PlotDataset::~PlotDataset() = default;

void plots::PlotDataset::initial_plot(const Dataset& data) {
    plots::PlotOptions options = data.plot_options;
    options.use_existing_axes = false;

    draw(data, options, canvas);
}

void plots::PlotDataset::plot(const Dataset& data) {
    plots::PlotOptions options = data.plot_options;
    options.use_existing_axes = true;

    draw(data, options, canvas);
}

void plots::PlotDataset::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotDataset::quick_plot(const Dataset& data, std::string path) {
    plots::PlotDataset plot(data);
    plot.save(path);
}

void plots::PlotDataset::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
}