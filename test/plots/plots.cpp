#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist//intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <fitter/HydrationFitter.h>
#include <settings/GeneralSettings.h>
#include <plots/All.h>

TEST_CASE("plot_dataset", "[manual]") {
    SimpleDataset data("test/files/2epe.dat");
    plots::PlotDataset::quick_plot(data, {}, settings::general::output + "test/plot/dataset.png");
}

TEST_CASE("plot_distance", "[manual]") {
    data::Molecule protein("test/files/2epe.pdb");
    protein.generate_new_hydration();
    auto data = protein.get_histogram();
    plots::PlotDistance::quick_plot(data.get(), settings::general::output + "test/plot/distance.png");
}

TEST_CASE("plot_histogram", "[manual]") {
    
}

TEST_CASE("plot_image", "[manual]") {
    
}

TEST_CASE("plot_intensity", "[manual]") {
    SimpleDataset data("test/files/2epe.dat");
    plots::PlotIntensity::quick_plot(data, {}, settings::general::output + "test/plot/dataset.png");    
}

TEST_CASE("plot_intensityfit", "[manual]") {
    data::Molecule protein("test/files/2epe.pdb");
    protein.generate_new_hydration();

    auto h = protein.get_histogram();
    fitter::HydrationFitter fitter("test/files/2epe.dat", std::move(h));
    auto result = fitter.fit();
    plots::PlotIntensityFit::quick_plot(result.get(), settings::general::output + "test/plot/intensityfit.png");
    plots::PlotIntensityFitResiduals::quick_plot(result.get(), settings::general::output + "test/plot/intensityfitresiduals.png");
}

TEST_CASE("plot_resolution", "[manual]") {
    
}