#include <plots/PlotDistance.h>
#include <histogram/DebyeLookupTable.h>
#include <utility/Settings.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::string;

plots::PlotDistance::~PlotDistance() = default;

plots::PlotDistance::PlotDistance(const hist::ScatteringHistogram& d, std::string path) {
    plot(d);
    save(path);
}

void plots::PlotDistance::plot(const hist::ScatteringHistogram& d) {
    PlotOptions options;
    options.xlabel = "Distance [Å]";
    options.ylabel = "Count";

    SimpleDataset p(d.p.data, d.q);
    SimpleDataset pp(d.p_pp.p.data, d.q);
    SimpleDataset ph(d.p_hp.p.data, d.q);
    SimpleDataset hh(d.p_hh.p.data, d.q);

    p.set_plot_options({{"color", "k"}, {"label", "total"}});
    pp.set_plot_options({{"color", "tab:orange"}, {"label", "atom-atom"}});
    ph.set_plot_options({{"color", "tab:green"}, {"label", "atom-water"}});
    hh.set_plot_options({{"color", "tab:blue"}, {"label", "water-water"}});

    ss << "PlotDistance\np\n"
        << p.to_string()
        << "\n"
        << p.get_plot_options().to_string()
        << "\npp\n"
        << pp.to_string()
        << "\n"
        << pp.get_plot_options().to_string()
        << "\nph\n"
        << ph.to_string()
        << "\n"
        << ph.get_plot_options().to_string()
        << "\nhh\n"
        << hh.to_string()
        << "\n"
        << hh.get_plot_options().to_string()
        << std::endl;
}