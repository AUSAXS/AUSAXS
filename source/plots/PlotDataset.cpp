/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotDataset.h>
#include <plots/PlotOptions.h>
#include <dataset/Dataset.h>
#include <dataset/Multiset.h>

using namespace plots;

PlotDataset::PlotDataset(const Dataset& data, const PlotOptions& options) {
    plot(data, options);
}

PlotDataset::~PlotDataset() = default;

PlotDataset& PlotDataset::vline(double x, const PlotOptions& options) {
    ss << "PlotVline\n" 
        << "  " << x << "\n"
        << options.to_string() 
        << std::endl;
    return *this;
}

PlotDataset& PlotDataset::hline(double y, const PlotOptions& options) {
    ss << "PlotHline\n" 
        << "  " << y << "\n"
        << options.to_string() 
        << std::endl;
    return *this;
}

PlotDataset& PlotDataset::plot(const Dataset& data, const PlotOptions& options) {
    ss << "PlotDataset\n" 
        << data.to_string() 
        << "\n"
        << options.to_string() 
        << std::endl;
    return *this;
}

void PlotDataset::quick_plot(const Dataset& data, const PlotOptions& options, const io::File& path) {
    PlotDataset plot(data, options);
    plot.save(path);
}