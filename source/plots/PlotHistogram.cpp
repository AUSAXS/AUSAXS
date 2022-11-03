#include <plots/PlotHistogram.h>
#include <utility/Utility.h>

plots::PlotHistogram::PlotHistogram(const hist::Histogram& h) {
    plot(h);
}

plots::PlotHistogram::~PlotHistogram() = default;

void plots::PlotHistogram::plot(const hist::Histogram& hist) {
    PlotOptions options;
    SimpleDataset p(hist.p.data, hist.axis.as_vector());

    ss << "PlotHistogram\n"
        << p.to_string()
        << "\n"
        << p.get_plot_options().to_string()
        << std::endl;
}

void plots::PlotHistogram::quick_plot(const hist::Histogram& hist, std::string path) {
    plots::PlotHistogram plot(hist);
    plot.save(path);
}