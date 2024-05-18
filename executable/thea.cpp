
#include <plots/All.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <constants/Constants.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>

#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    settings::em::sample_frequency = 5;

    em::ImageStack stack(argv[1]);
    Axis axis = {10, 20, 10};

    std::vector<hist::ScatteringProfile> histograms;
    for (double cutoff = axis.max; cutoff >= axis.min; cutoff -= axis.width()) {
        std::cout << "cutoff: " << cutoff << std::endl;
        stack.save(cutoff, "output/thea/protein_" + std::to_string(cutoff) + ".pdb");
        // histograms.push_back(stack.get_histogram(cutoff)->debye_transform());
    }

    plots::PlotDataset plot;
    for (const auto& h : histograms) {
        static double stagger = 1;
        plot.plot(h.as_dataset(), plots::PlotOptions({{"title", "Cutoff Scan"}, {"xlabel", "Distance"}, {"ylabel", "Frequency"}, {"yrange", std::vector{1e-4, 1e10}}, {"stagger", stagger}}));
        stagger*=10;
    }
    plot.save("output/thea/cutoff_scan.png");
}