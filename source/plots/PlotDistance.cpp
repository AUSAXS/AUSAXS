#include <plots/PlotDistance.h>
#include <plots/PlotDataset.h>
#include <hist/DebyeLookupTable.h>
#include <utility/Settings.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::string;

plots::PlotDistance::~PlotDistance() = default;

plots::PlotDistance::PlotDistance(const hist::ScatteringHistogram& d, std::string path) {
    quick_plot(d, path);
}

void plots::PlotDistance::quick_plot(const hist::ScatteringHistogram& d, std::string path) {
    PlotOptions options;
    options.xlabel = "Distance [Å]";
    options.ylabel = "Count";

    auto distances = d.axis.as_vector();
    SimpleDataset p(distances, d.p.data);
    SimpleDataset pp(distances, d.p_pp.p.data);
    SimpleDataset ph(distances, d.p_hp.p.data);
    SimpleDataset hh(distances, d.p_hh.p.data);

    p.add_plot_options("lines", {{"color", style::color::black}, {"legend", "total"}});
    pp.add_plot_options("lines", {{"color", style::color::orange}, {"legend", "atom-atom"}});
    ph.add_plot_options("lines", {{"color", style::color::green}, {"legend", "atom-water"}});
    hh.add_plot_options("lines", {{"color", style::color::blue}, {"legend", "water-water"}});

    plots::PlotDataset plot;
    plot.plot(p);
    plot.plot(pp);
    plot.plot(ph);
    plot.plot(hh);
    plot.save(path);
}
