#include <plots/PlotIntensityFit.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

plots::PlotIntensityFit::PlotIntensityFit(fitter::LinearFitter& fitter) : Plot() {
    auto graphs = fitter.plot();
    plot(graphs);
}

plots::PlotIntensityFit::PlotIntensityFit(const Fit& fit) : Plot() {
    plot(fit.figures);
}

plots::PlotIntensityFit::PlotIntensityFit(const std::shared_ptr<Fit> fit) : Plot() {
    plot(fit->figures);
}

plots::PlotIntensityFit::~PlotIntensityFit() = default;

void plots::PlotIntensityFit::quick_plot(const std::shared_ptr<Fit> fit, std::string path) {
    PlotIntensityFit plot(fit);
    plot.save(path);
}

void plots::PlotIntensityFit::plot(const Fit::Plots& graphs) {
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