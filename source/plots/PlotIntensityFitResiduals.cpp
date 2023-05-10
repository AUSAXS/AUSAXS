#include <plots/PlotIntensityFitResiduals.h>
#include <fitter/Fit.h>
#include <hist/ScatteringHistogram.h>
#include <fitter/LinearFitter.h>

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(fitter::LinearFitter& fitter) : Plot() {
    SimpleDataset graph = fitter.plot_residuals();
    plot(graph);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const fitter::Fit& fit) : Plot() {
    plot(fit.residuals);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const std::shared_ptr<fitter::Fit> fit) : Plot() {
    plot(fit->residuals);
}

plots::PlotIntensityFitResiduals::~PlotIntensityFitResiduals() = default;

void plots::PlotIntensityFitResiduals::quick_plot(const std::shared_ptr<fitter::Fit> fit, const io::File& path) {
    PlotIntensityFitResiduals plot(fit);
    plot.save(path);
}

void plots::PlotIntensityFitResiduals::plot(const SimpleDataset& graph) {
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