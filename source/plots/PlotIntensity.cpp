#include <fitter/Fit.h>
#include <plots/PlotIntensity.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

plots::PlotIntensity::PlotIntensity(const hist::ScatteringHistogram& d, std::string color = color::black) {
    plot(d, color);
}

plots::PlotIntensity::PlotIntensity(const SimpleDataset& d, std::string color = color::black) {
    plot(d, color);
}

plots::PlotIntensity::~PlotIntensity() = default;

void plots::PlotIntensity::plot(const SimpleDataset& data, std::string color = color::black) {
    plots::PlotOptions options = data.get_plot_options();
    options.xlabel = "q";
    options.ylabel = "Intensity";
    options.ylimits = data.span_y_positive();
    options.color = color;

    ss << "PlotIntensity\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
}

void plots::PlotIntensity::plot(const SimpleDataset& data, std::string color = color::black) {
    PlotOptions options(data.get_plot_options());
    options.use_existing_axes = true;
    options.color = color;

    ss << "PlotIntensity\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
}

void plots::PlotIntensity::plot(const std::shared_ptr<Fit> fit, const PlotOptions& plot_options) {
    plot(fit->figures.intensity);
}

void plots::PlotIntensity::plot_guinier_approx(const hist::ScatteringHistogram& data) {
    PlotOptions options = data.get_plot_options();
    options.legend = "Guinier approx";
    options.xlabel = "q [Å]";
    options.ylabel = "I";

    auto guinier = data.plot_guinier_approx();
    ss << "Guinier approx\n"
       << guinier.to_string()
       << "\n"
       << options.to_string()
       << std::endl;

    ss << "Gyration ratio\n"
        << sqrt(data.calc_guinier_gyration_ratio_squared())
        << std::endl;
}

void plots::PlotIntensity::quick_plot(const SimpleDataset& d, std::string path) {
    PlotIntensity plot(d);
    plot.save(path);
}

void plots::PlotIntensity::quick_plot(const hist::ScatteringHistogram& h, std::string path) {
    PlotIntensity plot(h);
    plot.save(path);
}