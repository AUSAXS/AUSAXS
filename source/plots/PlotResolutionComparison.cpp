#include <plots/PlotResolutionComparison.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::string, std::vector;

plots::PlotResolutionComparison::PlotResolutionComparison(Multiset data) {
    if (data.empty()) {throw except::size_error("Error in PlotResolutionComparison::PlotResolutionComparison: The given Multiset is empty!");}

    // change colors & plot raw figure before we scale the y-values
    data[0].add_plot_options("line", {{"logx", true}, {"logy", true}, {"xlabel", "q"}, {"ylabel", "Intensity"}, {"ylimit", Limit(1e-4, inf)}});
    plots::PlotDataset::quick_plot(data, "resolution_raw");

    // scale the y-values & plot
    data[0].scale_y(std::pow(1.3, data.size()));
    for (unsigned int i = 1; i < data.size(); i++) {
        data[i].scale_y(std::pow(1.3, data.size()-i));   // scale the y-axis so the plots becomes staggered
    }
    plots::PlotDataset::quick_plot(data, "resolution_staggered");
}

plots::PlotResolutionComparison::~PlotResolutionComparison() = default;

void plots::PlotResolutionComparison::quick_plot(const Multiset& data, std::string path) {
    plots::PlotResolutionComparison plot(data);
    plot.save(path);
}