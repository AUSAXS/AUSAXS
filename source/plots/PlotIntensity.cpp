#include <plots/PlotIntensity.h>
#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotIntensity::PlotIntensity(const hist::ScatteringProfile& d, style::Color color) {
    plot(d, color);
}

PlotIntensity::PlotIntensity(const SimpleDataset& d, style::Color color) {
    plot(d, color);
}

PlotIntensity::~PlotIntensity() = default;

PlotIntensity& PlotIntensity::plot(const SimpleDataset& data, style::Color color) {
    PlotOptions options = data.get_plot_options();
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    options.ylimits = data.span_y_positive();
    options.ylimits.max *= 1.1;
    options.logx = true;
    options.logy = true;
    options.color = color;

    ss << "PlotDataset\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
    return *this;
}

PlotIntensity& PlotIntensity::plot(const hist::ScatteringProfile& data, style::Color color) {
    PlotOptions options;
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
    return *this;
}

PlotIntensity& PlotIntensity::plot(const view_ptr<fitter::Fit> fit, style::Color color) {
    return plot(fit->figures.intensity, color);
}

PlotIntensity& PlotIntensity::plot_guinier_approx(const view_ptr<hist::CompositeDistanceHistogram> data) {
    PlotOptions options;
    options.legend = "Guinier approx";
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    options.logx = true;
    options.logy = true;

    throw except::not_implemented("PlotIntensity::plot_guinier_approx: Not implemented yet!");

    // auto guinier = data->plot_guinier_approx();
    // ss << "Guinier\n"
    //    << guinier.to_string()
    //    << "\n"
    //    << options.to_string()
    //    << std::endl;

    // ss << "Gyration_ratio\n"
    //     << sqrt(data->calc_guinier_gyration_ratio_squared())
    //     << std::endl;
    // return *this;
}

void PlotIntensity::quick_plot(const SimpleDataset& d, const io::File& path) {
    PlotIntensity plot(d);
    plot.save(path);
}

void PlotIntensity::quick_plot(const hist::ScatteringProfile& h, const io::File& path) {
    PlotIntensity plot(h);
    plot.save(path);
}