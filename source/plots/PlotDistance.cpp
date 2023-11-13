#include <plots/PlotDistance.h>
#include <plots/PlotDataset.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotDistance::~PlotDistance() = default;

PlotDistance::PlotDistance(const view_ptr<hist::ICompositeDistanceHistogram> d, const io::File& path) {
    quick_plot(d, path);
}

void PlotDistance::quick_plot(const view_ptr<hist::ICompositeDistanceHistogram> d, const io::File& path) {
    auto distances = d->get_axis().as_vector();
    SimpleDataset p(distances, d->get_counts());
    SimpleDataset pp(distances, d->get_aa_counts());
    SimpleDataset ph(distances, d->get_aw_counts());
    SimpleDataset hh(distances, d->get_ww_counts());

    PlotDataset plot;
    plot.plot(p,  plots::PlotOptions("lines", {{"color", style::color::black}, {"legend", "total"}, {"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}}));
    plot.plot(pp, plots::PlotOptions("lines", {{"color", style::color::orange}, {"legend", "atom-atom"}}));
    plot.plot(ph, plots::PlotOptions("lines", {{"color", style::color::green}, {"legend", "atom-water"}}));
    plot.plot(hh, plots::PlotOptions("lines", {{"color", style::color::blue}, {"legend", "water-water"}}));
    plot.save(path);
}
