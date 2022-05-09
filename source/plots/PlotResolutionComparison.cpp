#include <plots/PlotResolutionComparison.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TCanvas.h>

using std::string, std::vector;

plots::PlotResolutionComparison::PlotResolutionComparison(const Multiset& data, int color) : Plot() {
    gStyle->SetPalette(color);
    auto cols = TColor::GetPalette();
    unsigned int color_step = (cols.GetSize()-1)/data.size();
    prepare_canvas();

    Dataset data_copy = data[0];
    data_copy.set_plot_options(cols.At(0*color_step));
    initial_plot(data_copy);
    for (unsigned int i = 1; i < data.size(); i++) {
        data_copy = data[i];
        data_copy.scale_y(std::pow(1.3, i));               // scale the y-axis so the plots becomes staggered
        data_copy.set_plot_options(cols.At(i*color_step)); // change the color to be a nice gradient
        plot(data_copy);
    }
}

plots::PlotResolutionComparison::~PlotResolutionComparison() = default;

void plots::PlotResolutionComparison::initial_plot(Dataset& data) {
    auto cols = TColor::GetPalette();
    draw(data, canvas);
}

void plots::PlotResolutionComparison::plot(Dataset& data) {
    plots::PlotOptions options = data.plot_options;
    options.use_existing_axes = true;

    draw(data, options, canvas);
}

void plots::PlotResolutionComparison::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

void plots::PlotResolutionComparison::quick_plot(const Multiset& data, std::string path) {
    plots::PlotResolutionComparison plot(data);
    plot.save(path);
}

void plots::PlotResolutionComparison::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
}