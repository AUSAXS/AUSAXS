#include <plots/PlotDataset.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <memory.h>
#include <string.h>
#include <vector>

#include <TCanvas.h>

using std::string, std::vector;

template<typename T>
plots::PlotDataset::PlotDataset(const T& data) : Plot() {
    prepare_canvas();
    initial_plot(data);
}

plots::PlotDataset::PlotDataset(const Multiset& data) {
    if (data.empty()) {throw except::size_error("Error in PlotDataset::PlotDataset: The given Multiset is empty!");}

    prepare_canvas();
    draw(data, canvas);
}

plots::PlotDataset::~PlotDataset() = default;

template<typename T>
void plots::PlotDataset::initial_plot(const T& data) {
    plots::PlotOptions options = data.get_plot_options();
    options.use_existing_axes = false;

    draw(data, options, canvas);
}

template<typename T>
void plots::PlotDataset::plot(const T& data) {
    plots::PlotOptions options = data.get_plot_options();
    options.use_existing_axes = true;

    draw(data, options, canvas);
}

void plots::PlotDataset::save(std::string path) const {
    utility::create_directories(path);
    canvas->SaveAs(path.c_str());
}

template<typename T>
void plots::PlotDataset::quick_plot(const T& data, std::string path) {
    plots::PlotDataset plot(data);
    plot.save(path);
}

template void plots::PlotDataset::quick_plot(const Dataset& data, std::string path);
template void plots::PlotDataset::quick_plot(const SimpleDataset& data, std::string path);
template plots::PlotDataset::PlotDataset(const Dataset& data);
template plots::PlotDataset::PlotDataset(const SimpleDataset& data);
template void plots::PlotDataset::initial_plot(const Dataset& data);
template void plots::PlotDataset::initial_plot(const SimpleDataset& data);
template void plots::PlotDataset::plot(const Dataset& data);
template void plots::PlotDataset::plot(const SimpleDataset& data);