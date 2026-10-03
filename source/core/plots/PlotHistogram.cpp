// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <plots/PlotHistogram.h>

#include <dataset/Dataset.h>
#include <hist/Histogram.h>
#include <hist/ScatteringProfile.h>

using namespace ausaxs::plots;

PlotHistogram::PlotHistogram() = default;

PlotHistogram::~PlotHistogram() = default;

PlotHistogram::PlotHistogram(const hist::Histogram& h, const PlotOptions& options) {
    plot(h, options);
}

PlotHistogram& PlotHistogram::plot(const hist::Histogram& hist, const PlotOptions& options) {
    ss << "PlotHistogram\n"
        << hist.as_dataset().to_string()
        << "\n"
        << options.to_string()
        << std::endl;
    return *this;
}

PlotHistogram& PlotHistogram::plot(const hist::ScatteringProfile& profile, const PlotOptions& options) {
    ss << "PlotHistogram\n"
        << profile.as_dataset().to_string()
        << "\n"
        << options.to_string()
        << std::endl;
    return *this;
}

void PlotHistogram::quick_plot(const hist::Histogram& hist, const PlotOptions& options, const io::File& path) {
    PlotHistogram plot(hist, options);
    plot.save(path);
}