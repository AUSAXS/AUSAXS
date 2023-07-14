#include <plots/PlotDataset.h>
#include <plots/PlotOptions.h>
#include <dataset/Dataset.h>
#include <dataset/Multiset.h>

using namespace plots;

template<plots::DatasetType T>
PlotDataset::PlotDataset(const T& data) {
    plot(data);
}

PlotDataset::PlotDataset(const Multiset& data) {
    if (data.empty()) {throw except::size_error("PlotDataset::PlotDataset: The given Multiset is empty!");}

    for (const auto& d : data) {
        plot(d);
    }
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

template<plots::DatasetType T>
PlotDataset& PlotDataset::plot(const T& data) {
    ss << "PlotDataset\n" 
        << data.to_string() 
        << "\n"
        << data.get_plot_options().to_string() 
        << std::endl;
    return *this;
}

template<plots::DatasetType T>
void PlotDataset::quick_plot(const T& data, const io::File& path) {
    PlotDataset plot(data);
    plot.save(path);
}

template void PlotDataset::quick_plot(const Dataset2D& data, const io::File& path);
template void PlotDataset::quick_plot(const SimpleDataset& data, const io::File& path);
template PlotDataset::PlotDataset(const Dataset2D& data);
template PlotDataset::PlotDataset(const SimpleDataset& data);
template PlotDataset& PlotDataset::plot(const Dataset2D& data);
template PlotDataset& PlotDataset::plot(const SimpleDataset& data);