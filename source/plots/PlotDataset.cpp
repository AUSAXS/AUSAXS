#include <plots/PlotDataset.h>
#include <Exceptions.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

plots::PlotDataset::PlotDataset(const Dataset& data) : Plot() {
    prepare_canvas();
    initial_plot(data);
}

void plots::PlotDataset::initial_plot(const Dataset& data) {
    std::shared_ptr<TGraph> graph = data.plot();
    plots::PlotOptions options = data.plot_options;
    options.xlabel = "cutoff";
    options.ylabel = "chi2";

    draw(graph, options);
}

void plots::PlotDataset::plot(const Dataset& data) {
    std::shared_ptr<TGraph> graph = data.plot();
    PlotOptions options(data.plot_options);
    options.use_existing_axes = true;

    draw(graph, options);
}

void plots::PlotDataset::save(std::string path) const {
    canvas->SaveAs(path.c_str());
}

void plots::PlotDataset::prepare_canvas() {
    canvas = std::make_unique<TCanvas>("PlotDatasetCanvas", "canvas", 600, 600);
}