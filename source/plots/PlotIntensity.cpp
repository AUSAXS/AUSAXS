#include <plots/PlotIntensity.h>
#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <dataset/SimpleDataset.h>

using namespace plots;

PlotIntensity::PlotIntensity(const hist::ScatteringProfile& d, const PlotOptions& options) {
    plot(d, options);
}

PlotIntensity::PlotIntensity(const SimpleDataset& d, const PlotOptions& options) {
    plot(d, options);
}

PlotIntensity::~PlotIntensity() = default;

PlotIntensity& PlotIntensity::plot(const SimpleDataset& data, const PlotOptions& temp) {
    PlotOptions options = temp;
    options.xlabel = "$q$ [$\\AA^{-1}$]";
    options.ylabel = "$I$ [arb]";
    if (options.ylimits.empty()) {
        options.ylimits = data.span_y_positive();
        options.ylimits.max *= 1.1;
    }
    options.logx = true;
    options.logy = true;

    ss << "PlotDataset\n"
       << data.to_string()
       << "\n"
       << options.to_string()
       << std::endl;
    return *this;
}

PlotIntensity& PlotIntensity::plot(const hist::ScatteringProfile& data, const PlotOptions& temp) {
    PlotOptions options = temp;
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

PlotIntensity& PlotIntensity::plot(std::observer_ptr<fitter::Fit> fit, const PlotOptions& options) {
    return plot(fit->figures.intensity, options);
}

PlotIntensity& PlotIntensity::plot_guinier_approx(std::observer_ptr<hist::ICompositeDistanceHistogram> data) {
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

void PlotIntensity::quick_plot(const SimpleDataset& d, const PlotOptions& options, const io::File& path) {
    PlotIntensity plot(d, options);
    plot.save(path);
}

void PlotIntensity::quick_plot(const hist::ScatteringProfile& h, const PlotOptions& options, const io::File& path) {
    PlotIntensity plot(h, options);
    plot.save(path);
}