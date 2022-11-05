#include <plots/PlotIntensityFitResiduals.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(SimpleIntensityFitter& fitter) : Plot() {
    SimpleDataset graph = fitter.plot_residuals();
    plot(graph);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const Fit& fit) : Plot() {
    plot(fit.residuals);
}

plots::PlotIntensityFitResiduals::PlotIntensityFitResiduals(const std::shared_ptr<Fit> fit) : Plot() {
    plot(fit->residuals);
}

plots::PlotIntensityFitResiduals::~PlotIntensityFitResiduals() = default;

void plots::PlotIntensityFitResiduals::quick_plot(const std::shared_ptr<Fit> fit, std::string path) {
    PlotIntensityFitResiduals plot(fit);
    plot.save(path);
}

void plots::PlotIntensityFitResiduals::plot(const SimpleDataset graph) {
    PlotOptions options("markers", {{"color", style::color::orange}, {"title", "Residuals"}, {"xlabel", "q"}, {"ylabel", "Residual"}});
    PlotOptions line("line", {{"color", style::color::black}});

    ss << "PlotIntensityFitResiduals\n"
        << graph.to_string()
        << "\n"
        << options.to_string()
        << std::endl
        << "hline\n"
        << 0
        << "\n"
        << line.to_string()
        << std::endl;
}