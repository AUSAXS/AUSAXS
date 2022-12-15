#include <fitter/Fit.h>
#include <plots/PlotIntensity.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <memory.h>
#include <string.h>
#include <vector>

plots::PlotIntensity::PlotIntensity(const hist::ScatteringHistogram& d, style::Color color) {
    plot(d, color);
}

plots::PlotIntensity::PlotIntensity(const SimpleDataset& d, style::Color color) {
    plot(d, color);
}

plots::PlotIntensity::~PlotIntensity() = default;

void plots::PlotIntensity::plot(const SimpleDataset& data, style::Color color) {
    plots::PlotOptions options = data.get_plot_options();
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    options.ylimits = data.span_y_positive();
    options.logx = true;
    options.logy = true;
    options.color = color;

    ss << "PlotDataset\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
}

void plots::PlotIntensity::plot(const hist::ScatteringHistogram& data, style::Color color) {
    PlotOptions options(data.get_plot_options());
    options.color = color;
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    options.logx = true;
    options.logy = true;

    ss << "PlotDataset\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
}

void plots::PlotIntensity::plot(const std::shared_ptr<Fit> fit, style::Color color) {
    plot(fit->figures.intensity, color);
}

void plots::PlotIntensity::plot_guinier_approx(const hist::ScatteringHistogram& data) {
    PlotOptions options = data.get_plot_options();
    options.legend = "Guinier approx";
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    options.logx = true;
    options.logy = true;

    auto guinier = data.plot_guinier_approx();
    ss << "Guinier\n"
       << guinier.to_string()
       << "\n"
       << options.to_string()
       << std::endl;

    ss << "Gyration_ratio\n"
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