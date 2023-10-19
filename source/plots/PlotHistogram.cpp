#include <plots/PlotHistogram.h>
#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotHistogram::PlotHistogram() = default;

PlotHistogram::~PlotHistogram() = default;

PlotHistogram::PlotHistogram(const hist::Histogram& h) {
    plot(h);
}

PlotHistogram& PlotHistogram::plot(const hist::Histogram& hist, const plots::PlotOptions& options) {
    SimpleDataset p(hist.get_axis().as_vector(), hist.get_counts());

    ss << "PlotHistogram\n"
        << p.to_string()
        << "\n"
        << options.to_string()
        << std::endl;
    return *this;
}

PlotHistogram& PlotHistogram::plot(const hist::Histogram& hist) {
    return plot(hist, hist.get_plot_options());
}

void PlotHistogram::quick_plot(const hist::Histogram& hist, const io::File& path) {
    PlotHistogram plot(hist);
    plot.save(path);
}