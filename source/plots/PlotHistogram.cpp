#include <plots/PlotHistogram.h>
#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

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