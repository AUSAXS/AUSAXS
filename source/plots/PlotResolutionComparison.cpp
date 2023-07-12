#include <plots/PlotResolutionComparison.h>
#include <plots/PlotDataset.h>
#include <dataset/Multiset.h>
#include <utility/Exceptions.h>

using namespace plots;

PlotResolutionComparison::PlotResolutionComparison(Multiset& data) {
    if (data.empty()) {throw except::size_error("PlotResolutionComparison::PlotResolutionComparison: The given Multiset is empty!");}

    // change colors & plot raw figure before we scale the y-values
    data[0].add_plot_options("line", {{"logx", true}, {"logy", true}, {"xlabel", "q"}, {"ylabel", "Intensity"}, {"ylimit", Limit(1e-4, inf)}});
    PlotDataset::quick_plot(data, "resolution_raw");

    // scale the y-values & plot
    data[0].scale_y(std::pow(1.3, data.size()));
    for (unsigned int i = 1; i < data.size(); i++) {
        data[i].scale_y(std::pow(1.3, data.size()-i));   // scale the y-axis so the plots becomes staggered
    }
    PlotDataset::quick_plot(data, "resolution_staggered");
}

PlotResolutionComparison::~PlotResolutionComparison() = default;

void PlotResolutionComparison::quick_plot(Multiset& data, const io::File& path) {
    PlotResolutionComparison plot(data);
    plot.save(path);
}