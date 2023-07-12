#include <plots/PlotHistogram.h>
#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotHistogram::PlotHistogram(const hist::Histogram& h) {
    plot(h);
}

PlotHistogram::~PlotHistogram() = default;

PlotHistogram& PlotHistogram::plot(const hist::Histogram& hist) {
    PlotOptions options;
    SimpleDataset p(hist.p.data, hist.axis.as_vector());

    ss << "PlotHistogram\n"
        << p.to_string()
        << "\n"
        << p.get_plot_options().to_string()
        << std::endl;
    return *this;
}

void PlotHistogram::quick_plot(const hist::Histogram& hist, const io::File& path) {
    PlotHistogram plot(hist);
    plot.save(path);
}