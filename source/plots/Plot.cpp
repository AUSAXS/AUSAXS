#include <plots/Plot.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <hist/Histogram.h>

#include <fstream>

void plots::Plot::save(std::string path) const {
    utility::create_directory(path);
    std::ofstream file(path + ".plot");
    if (!file.is_open()) {throw except::io_error("PlotDataset::quick_plot: Could not open file " + path + " for writing!");}
    file << ss.str();
    file.close();
}