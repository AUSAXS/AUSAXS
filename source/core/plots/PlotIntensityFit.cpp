/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotIntensityFit.h>
#include <plots/PlotOptions.h>
#include <fitter/LinearFitter.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

using namespace plots;

PlotIntensityFit::PlotIntensityFit(fitter::LinearFitter& fitter) : Plot() {
    auto graphs = fitter.plot();
    plot(graphs);
}

PlotIntensityFit::PlotIntensityFit(const fitter::FitResult& fit) : Plot() {
    plot(fit.info);
}

PlotIntensityFit::PlotIntensityFit(observer_ptr<fitter::FitResult> fit) : Plot() {
    plot(fit->info);
}

PlotIntensityFit::~PlotIntensityFit() = default;

void PlotIntensityFit::quick_plot(observer_ptr<fitter::FitResult> fit, const io::File& path) {
    PlotIntensityFit plot(fit);
    plot.save(path);
}

void PlotIntensityFit::plot(const fitter::FitResult::FitInfo& graphs) {
    PlotOptions options_data, options_interpolated, options_intensity;
    options_data.set("errors", {{"color", style::color::orange}, {"title", "Fit"}, {"xlabel", "$q$ [$\\AA^{-1}$]"}, {"ylabel", "$I$ [arb]"}, {"logy", true}, {"logx", true}});
    options_interpolated.set("markers", {{"color", style::color::black}});
    options_intensity.set("line", {{"color", style::color::black}});

    ss << "PlotDataset\n"
       << graphs.dataset.to_string()
       << "\n"
       << options_data.to_string()
       << "\nPlotDataset\n"
       << graphs.fitted_intensity_interpolated.to_string()
       << "\n"
       << options_interpolated.to_string()
       << "\nPlotDataset\n"
       << graphs.fitted_intensity.to_string()
       << "\n"
       << options_intensity.to_string()
       << std::endl;
}