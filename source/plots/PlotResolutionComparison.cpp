#include <plots/PlotResolutionComparison.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TCanvas.h>

using std::string, std::vector;

plots::PlotResolutionComparison::PlotResolutionComparison(Multiset data, int color) : Plot() {
    // note that the data is passed as a copy, so we're free to change the plot options
    gStyle->SetPalette(color);
    auto cols = TColor::GetPalette();
    unsigned int color_step = (cols.GetSize()-1)/data.size();
    prepare_canvas();

    // change colors & plot raw figure before we scale the y-values
    data[0].add_plot_options({{"logx", true}, {"logy", true}, {"xlabel", "q"}, {"ylabel", "Intensity"}});
    for (unsigned int i = 0; i < data.size(); i++) {
        data[i].add_plot_options(cols.At(i*color_step));
    }
    raw = std::make_unique<PlotDataset>(data);
    canvas->cd();

    // scale the y-values & plot
    data[0].scale_y(std::pow(1.3, data.size()));
    initial_plot(data[0]);
    for (unsigned int i = 1; i < data.size(); i++) {
        data[i].scale_y(std::pow(1.3, data.size()-i));   // scale the y-axis so the plots becomes staggered
        plot(data[i]);
    }
}

plots::PlotResolutionComparison::~PlotResolutionComparison() = default;

void plots::PlotResolutionComparison::initial_plot(Dataset& data) {
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

    raw->save(utility::stem_append(path, "_raw"));
}

void plots::PlotResolutionComparison::quick_plot(const Multiset& data, std::string path) {
    plots::PlotResolutionComparison plot(data);
    plot.save(path);
}

void plots::PlotResolutionComparison::prepare_canvas() {
    canvas = std::make_unique<TCanvas>(utility::uid("canvas").c_str(), "canvas", 600, 600);
}