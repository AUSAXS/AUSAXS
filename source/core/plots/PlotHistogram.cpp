/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotHistogram.h>
#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

using namespace ausaxs::plots;

PlotHistogram::PlotHistogram() = default;

PlotHistogram::~PlotHistogram() = default;

PlotHistogram::PlotHistogram(const hist::Histogram& h, const PlotOptions& options) {
    plot(h, options);
}

PlotHistogram& PlotHistogram::plot(const hist::Histogram& hist, const PlotOptions& options) {
    SimpleDataset p(hist.get_axis().as_vector(), hist.get_counts());

    ss << "PlotHistogram\n"
        << p.to_string()
        << "\n"
        << options.to_string()
        << std::endl;
    return *this;
}

void PlotHistogram::quick_plot(const hist::Histogram& hist, const PlotOptions& options, const io::File& path) {
    PlotHistogram plot(hist, options);
    plot.save(path);
}