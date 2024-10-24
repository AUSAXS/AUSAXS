/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotDistance.h>
#include <plots/PlotDataset.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/Distribution1D.h>
#include <dataset/SimpleDataset.h>

using namespace ausaxs::plots;

PlotDistance::~PlotDistance() = default;

PlotDistance::PlotDistance(observer_ptr<hist::DistanceHistogram> d, const io::File& path) {
    quick_plot(d, path);
}

void PlotDistance::quick_plot(observer_ptr<hist::DistanceHistogram> d, const io::File& path) {
    const auto& distances = d->get_d_axis();
    SimpleDataset p(distances, d->get_total_counts());

    PlotDataset plot;
    plot.plot(p,  plots::PlotOptions("lines", {{"color", style::color::black}, {"legend", "total"}, {"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}}));
    if (auto cast = dynamic_cast<hist::ICompositeDistanceHistogram*>(d)) {
        SimpleDataset pp(distances, cast->get_aa_counts());
        SimpleDataset ph(distances, cast->get_aw_counts());
        SimpleDataset hh(distances, cast->get_ww_counts());

        plot.plot(pp, plots::PlotOptions("lines", {{"color", style::color::orange}, {"legend", "atom-atom"}}));
        plot.plot(ph, plots::PlotOptions("lines", {{"color", style::color::green}, {"legend", "atom-water"}}));
        plot.plot(hh, plots::PlotOptions("lines", {{"color", style::color::blue}, {"legend", "water-water"}}));
    }
    plot.save(path);
}