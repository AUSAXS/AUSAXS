#include <plots/PlotDataset.h>

template<plots::DatasetType T>
plots::PlotDataset::PlotDataset(const T& data) {
    plot(data);
}

plots::PlotDataset::PlotDataset(const Multiset& data) {
    if (data.empty()) {throw except::size_error("PlotDataset::PlotDataset: The given Multiset is empty!");}

    for (const auto& d : data) {
        plot(d);
    }
}

plots::PlotDataset::~PlotDataset() = default;

template<plots::DatasetType T>
void plots::PlotDataset::plot(const T& data) {
    ss << "PlotDataset\n" 
        << data.to_string() 
        << "\n"
        << data.get_plot_options().to_string() 
        << std::endl;
}

template<plots::DatasetType T>
void plots::PlotDataset::quick_plot(const T& data, const io::File& path) {
    PlotDataset plot(data);
    plot.save(path);
}

template void plots::PlotDataset::quick_plot(const Dataset2D& data, const io::File& path);
template void plots::PlotDataset::quick_plot(const SimpleDataset& data, const io::File& path);
template plots::PlotDataset::PlotDataset(const Dataset2D& data);
template plots::PlotDataset::PlotDataset(const SimpleDataset& data);
template void plots::PlotDataset::plot(const Dataset2D& data);
template void plots::PlotDataset::plot(const SimpleDataset& data);