#include <plots/PlotDistance.h>
#include <plots/PlotDataset.h>
#include <hist/CompositeDistanceHistogram.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotDistance::~PlotDistance() = default;

PlotDistance::PlotDistance(const hist::CompositeDistanceHistogram* const d, const io::File& path) {
    quick_plot(d, path);
}

void PlotDistance::quick_plot(const hist::CompositeDistanceHistogram* const d, const io::File& path) {
    auto distances = d->get_axis().as_vector();
    SimpleDataset p(distances, d->p.data);
    SimpleDataset pp(distances, d->get_pp_histogram());
    SimpleDataset ph(distances, d->get_hp_histogram());
    SimpleDataset hh(distances, d->get_hh_histogram());

    p.add_plot_options("lines", {{"color", style::color::black}, {"legend", "total"}, {"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
    pp.add_plot_options("lines", {{"color", style::color::orange}, {"legend", "atom-atom"}});
    ph.add_plot_options("lines", {{"color", style::color::green}, {"legend", "atom-water"}});
    hh.add_plot_options("lines", {{"color", style::color::blue}, {"legend", "water-water"}});

    PlotDataset plot;
    plot.plot(p);
    plot.plot(pp);
    plot.plot(ph);
    plot.plot(hh);
    plot.save(path);
}
