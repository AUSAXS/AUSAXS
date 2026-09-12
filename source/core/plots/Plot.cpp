// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <plots/Plot.h>

#include <io/File.h>
#include <utility/Exceptions.h>

#include <fstream>

using namespace ausaxs;

void plots::Plot::save(const io::File& path) const {
    path.directory().create();
    std::ofstream output(path.str() + ".plot");
    if (!output.is_open()) {throw except::io_error("PlotDataset::quick_plot: Could not open file " + path.str() + " for writing!");}
    output << ss.str();
}