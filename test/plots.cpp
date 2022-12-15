#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <plots/all.h>
#include <fitter/IntensityFitter.h>

TEST_CASE("plot_dataset", "[plots],[manual]") {
    SimpleDataset data("test/files/2epe.dat");
    plots::PlotDataset::quick_plot(data, "temp/plot/dataset.png");
}

TEST_CASE("plot_distance", "[plots],[manual]") {
    Protein protein("test/files/2epe.pdb");
    protein.generate_new_hydration();
    auto data = protein.get_histogram();
    plots::PlotDistance::quick_plot(data, "temp/plot/distance.png");
}

TEST_CASE("plot_histogram", "[plots],[manual]") {
    
}

TEST_CASE("plot_image", "[plots],[manual]") {
    
}

TEST_CASE("plot_intensity", "[plots],[manual]") {
    SimpleDataset data("test/files/2epe.dat");
    plots::PlotIntensity::quick_plot(data, "temp/plot/dataset.png");    
}

TEST_CASE("plot_intensityfit", "[plots],[manual]") {
    Protein protein("test/files/2epe.pdb");
    protein.generate_new_hydration();

    hist::ScatteringHistogram h = protein.get_histogram();
    IntensityFitter fitter("test/files/2epe.dat", h);
    std::shared_ptr<Fit> result = fitter.fit();
    plots::PlotIntensityFit::quick_plot(result, "temp/plot/intensityfit.png");
    plots::PlotIntensityFitResiduals::quick_plot(result, "temp/plot/intensityfitresiduals.png");
}

TEST_CASE("plot_resolution", "[plots],[manual]") {
    
}