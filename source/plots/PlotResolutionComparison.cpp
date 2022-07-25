#include <plots/PlotResolutionComparison.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TCanvas.h>

using std::string, std::vector;

plots::PlotResolutionComparison::PlotResolutionComparison(Multiset data, int color) : Plot() {
    if (data.empty()) {throw except::size_error("Error in PlotResolutionComparison::PlotResolutionComparison: The given Multiset is empty!");}

    // note that the data is passed as a copy, so we're free to change the plot options
    gStyle->SetPalette(color);
    auto cols = TColor::GetPalette();
    unsigned int color_step = (cols.GetSize()-1)/data.size();

    // change colors & plot raw figure before we scale the y-values
    data[0].add_plot_options({{"logx", true}, {"logy", true}, {"xlabel", "q"}, {"ylabel", "Intensity"}, {"ylimit", Limit(1e-4, inf)}});
    for (unsigned int i = 0; i < data.size(); i++) {
        data[i].set_plot_color(cols.At(i*color_step));
    }
    raw = std::make_unique<PlotDataset>(data);

    // scale the y-values & plot
    data[0].scale_y(std::pow(1.3, data.size()));
    for (unsigned int i = 1; i < data.size(); i++) {
        data[i].scale_y(std::pow(1.3, data.size()-i));   // scale the y-axis so the plots becomes staggered
    }
    staggered = std::make_unique<PlotDataset>(data);
}

plots::PlotResolutionComparison::~PlotResolutionComparison() = default;

void plots::PlotResolutionComparison::save(std::string path) const {
    utility::create_directories(path);
    staggered->save(path.c_str());
    raw->save(utility::stem_append(path, "_raw"));
}

void plots::PlotResolutionComparison::quick_plot(const Multiset& data, std::string path) {
    plots::PlotResolutionComparison plot(data);
    plot.save(path);
}