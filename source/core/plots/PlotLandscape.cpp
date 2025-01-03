/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotLandscape.h>
#include <mini/detail/Landscape.h>
#include <mini/detail/Evaluation.h>

using namespace ausaxs::plots;

PlotLandscape::PlotLandscape(const mini::Landscape& data, const PlotOptions& options, const io::File& path) {
    ss << "PlotLandscape\n" 
        << data.to_string() 
        << "\n"
        << options.to_string() 
        << std::endl;

    save(path);
}

PlotLandscape::~PlotLandscape() = default;

void PlotLandscape::quick_plot(const mini::Landscape& data, const PlotOptions& options, const io::File& path) {
    PlotLandscape plot(data, options, path);
}