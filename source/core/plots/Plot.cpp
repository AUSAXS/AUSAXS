/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/Plot.h>
#include <utility/Exceptions.h>
#include <hist/Histogram.h>
#include <io/File.h>

#include <fstream>

using namespace ausaxs;

void plots::Plot::save(const io::File& path) const {
    path.directory().create();
    std::ofstream output(path + ".plot");
    if (!output.is_open()) {throw except::io_error("PlotDataset::quick_plot: Could not open file " + path + " for writing!");}
    output << ss.str();
}