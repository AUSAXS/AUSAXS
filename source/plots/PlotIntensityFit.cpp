#include <plots/PlotIntensityFit.h>
#include <plots/PlotOptions.h>
#include <fitter/Fit.h>
#include <fitter/LinearFitter.h>
#include <fitter/FitPlots.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

using namespace plots;

PlotIntensityFit::PlotIntensityFit(fitter::LinearFitter& fitter) : Plot() {
    auto graphs = fitter.plot();
    plot(graphs);
}

PlotIntensityFit::PlotIntensityFit(const fitter::Fit& fit) : Plot() {
    plot(fit.figures);
}

PlotIntensityFit::PlotIntensityFit(const std::shared_ptr<fitter::Fit> fit) : Plot() {
    plot(fit->figures);
}

PlotIntensityFit::~PlotIntensityFit() = default;

void PlotIntensityFit::quick_plot(const std::shared_ptr<fitter::Fit> fit, const io::File& path) {
    PlotIntensityFit plot(fit);
    plot.save(path);
}

void PlotIntensityFit::plot(const fitter::FitPlots& graphs) {
    PlotOptions options_data, options_interpolated, options_intensity;
    options_data.set("errors", {{"color", style::color::orange}, {"title", "Fit"}, {"xlabel", "$q$ [$\\AA^{-1}$]"}, {"ylabel", "$I$ [arb]"}, {"logy", true}, {"logx", true}});
    options_interpolated.set("markers", {{"color", style::color::black}});
    options_intensity.set("line", {{"color", style::color::black}});

    ss << "PlotDataset\n"
       << graphs.data.to_string()
       << "\n"
       << options_data.to_string()
       << "\nPlotDataset\n"
       << graphs.intensity_interpolated.to_string()
       << "\n"
       << options_interpolated.to_string()
       << "\nPlotDataset\n"
       << graphs.intensity.to_string()
       << "\n"
       << options_intensity.to_string()
       << std::endl;
}