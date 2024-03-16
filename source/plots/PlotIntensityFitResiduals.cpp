/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/PlotIntensityFitResiduals.h>
#include <fitter/Fit.h>
#include <fitter/LinearFitter.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

using namespace plots;

PlotIntensityFitResiduals::PlotIntensityFitResiduals(fitter::LinearFitter& fitter) : Plot() {
    SimpleDataset graph = fitter.plot_residuals();
    plot(graph);
}

PlotIntensityFitResiduals::PlotIntensityFitResiduals(const fitter::Fit& fit) : Plot() {
    plot(fit.residuals);
}

PlotIntensityFitResiduals::PlotIntensityFitResiduals(observer_ptr<fitter::Fit> fit) : Plot() {
    plot(fit->residuals);
}

PlotIntensityFitResiduals::~PlotIntensityFitResiduals() = default;

void PlotIntensityFitResiduals::quick_plot(observer_ptr<fitter::Fit> fit, const io::File& path) {
    PlotIntensityFitResiduals plot(fit);
    plot.save(path);
}

void PlotIntensityFitResiduals::plot(const SimpleDataset& graph) {
    PlotOptions options("markers", {{"color", style::color::orange}, {"title", "Residuals"}, {"xlabel", "$q$ [$\\AA^{-1}$]"}, {"ylabel", "Residual"}, {"logx", true}});
    PlotOptions line("line", {{"color", style::color::black}});

    ss << "PlotDataset\n"
        << graph.to_string()
        << "\n"
        << options.to_string()
        << std::endl
        << "PlotHline\n"
        << 0
        << "\n\n"
        << line.to_string()
        << std::endl;
}