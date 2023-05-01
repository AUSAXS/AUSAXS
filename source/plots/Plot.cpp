#include <plots/Plot.h>
#include <utility/Exceptions.h>
#include <hist/Histogram.h>

#include <fstream>

void plots::Plot::save(const io::File& path) const {
    path.directory().create();
    std::ofstream output(path + ".plot");
    if (!output.is_open()) {throw except::io_error("PlotDataset::quick_plot: Could not open file " + path + " for writing!");}
    output << ss.str();
}