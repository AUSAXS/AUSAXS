#include <CLI/CLI.hpp>

#include <plots/All.h>
#include <em/ImageStack.h>
#include <data/Protein.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <fitter/LinearFitter.h>
#include <fitter/HydrationFitter.h>
#include <settings/All.h>
#include <utility/Constants.h>
#include <mini/dlibMinimizer.h>
#include <mini/detail/Parameter.h>

#include <iostream>
#include <unordered_set>
#include <tuple>

int main(int argc, char const *argv[]) {
    settings::protein::use_effective_charge = false;
    settings::em::mass_axis = true;
    settings::em::hydrate = false;
    settings::fit::verbose = true;

    io::ExistingFile mfile, mapfile1, mapfile2, mapfile3, settings;
    double cutoff1, cutoff2, cutoff3;
    CLI::App app{"Fit an EM map to a SAXS measurement."};
    app.add_option("input-exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    app.add_option("input-map1", mapfile1, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-map2", mapfile2, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-map3", mapfile3, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("cutoff1", cutoff1, "Cutoff for the first map.")->required();
    app.add_option("cutoff2", cutoff2, "Cutoff for the second map.")->required();
    app.add_option("cutoff3", cutoff3, "Cutoff for the third map.")->required();
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    // app.add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Whether to hydrate the protein before fitting.");
    // app.add_flag("--fixed-weight,!--no-fixed-weight", settings::em::fixed_weights, "Whether to use a fixed weight for the fit.");
    CLI11_PARSE(app, argc, argv);
    
    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }
    settings::general::output = "output/gregers/" + mfile.stem() + "/";
    std::filesystem::create_directories(settings::general::output);
    std::ofstream info(settings::general::output + "log.txt");
    info << "Fitting " << mfile << " with\n" 
            << "\t" << mapfile1.stem() << "\t" << cutoff1 << "\n"
            << "\t" << mapfile2.stem() << "\t" << cutoff2 << "\n"
            << "\t" << mapfile3.stem() << "\t" << cutoff3 << std::endl;
    info.close();

    //#####################################//
    //############## SAMPLER ##############//
    //#####################################//

    double area_factor = 0.2;
    unsigned int bins = 10;
    Dataset data(bins, 3);

    settings::em::fixed_weights = true;
    settings::em::sample_frequency = 2;
    em::ImageStack map1(mapfile1);
    std::vector<hist::ScatteringHistogram> histograms1;
    std::cout << "PREPARATION:\n\tPreparing map1" << std::endl;
    {
        Axis range(map1.from_level(cutoff1*(1 - area_factor)), map1.from_level(cutoff1*(1 + area_factor)), bins);
        for (unsigned int i = 0; i < bins; ++i) {
            double cutoff = range.min + i*range.width();
            auto h = map1.get_protein(cutoff)->get_histogram();
            h.extend_axis(1000);
            histograms1.push_back(h);
            data(i, 0) = map1.to_level(cutoff);
        }
    }

    settings::em::fixed_weights = false;
    settings::em::sample_frequency = 1;
    em::ImageStack map2(mapfile2);
    std::vector<hist::ScatteringHistogram> histograms2;
    std::cout << "\tPreparing map2" << std::endl;
    {
        Axis range(map2.from_level(cutoff2*(1 - area_factor)), map2.from_level(cutoff2*(1 + area_factor)), bins);
        for (unsigned int i = 0; i < bins; ++i) {
            double cutoff = range.min + i*range.width();
            auto h = map2.get_protein(cutoff)->get_histogram();
            h.extend_axis(1000);
            histograms2.push_back(h);
            data(i, 1) = map2.to_level(cutoff);
        }
    }

    settings::em::sample_frequency = 1;
    em::ImageStack map3(mapfile3);
    std::vector<hist::ScatteringHistogram> histograms3;
    std::cout << "\tPreparing map3" << std::endl;
    {
        Axis range(map3.from_level(cutoff3*(1 - area_factor)), map3.from_level(cutoff3*(1 + area_factor)), bins);
        for (unsigned int i = 0; i < bins; ++i) {
            double cutoff = range.min + i*range.width();
            auto h = map3.get_protein(cutoff)->get_histogram();
            h.extend_axis(1000);
            histograms3.push_back(h);
            data(i, 2) = map3.to_level(cutoff);
        }
    }
    std::cout << "Finished map preparation." << std::endl;

    //#####################################//
    //############## FITTER ##############//
    //#####################################//

    SimpleDataset dataset(mfile);
    std::vector<mini::Result> fits;
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> indices;
    fits.reserve(histograms1.size()*histograms2.size()*histograms3.size());
    indices.reserve(histograms1.size()*histograms2.size()*histograms3.size());
    for (unsigned int i = 0; i < histograms1.size(); ++i) {
        auto& h1 = histograms1[i];
        for (unsigned int j = 0; j < histograms2.size(); ++j) {
            auto& h2 = histograms2[j];
            for (unsigned int k = 0; k < histograms3.size(); ++k) {
                auto& h3 = histograms3[k];

                // 2 fit params: w1, w2
                auto func = [&] (std::vector<double> x) -> double {
                    double w1 = x[0];
                    double w2 = x[1]*(1-w1); // w2 = w2*(1-w1)
                    double w3 = 1 - w1 - w2; // w3 = (1-w1)*(1-w2)

                    hist::ScatteringHistogram composite = h1*w1 + h2*w2 + h3*w3;
                    return fitter::LinearFitter(dataset, composite).fit_chi2_only();
                };

                mini::dlibMinimizer<mini::type::BFGS> minimizer(func, std::vector{mini::Parameter{"x1", 0.33, {0, 1}}, mini::Parameter{"x2", 0.5, {0, 1}}});
                minimizer.set_max_evals(100);
                auto res = minimizer.minimize();
                std::cout << "(" << data.index(i, 0) << ", " << data.index(j, 1) << ", " << data.index(k, 2) << "): (x1, x2) = (" << res.get_parameter("x1").value << ", " << res.get_parameter("x2").value << ") with value "<< res.fval << std::endl;
                fits.push_back(res);
                indices.push_back(std::make_tuple(i, j, k));
            }
        }
    }

    //#####################################//
    //############## REPORT ##############//
    //#####################################//

    auto best = std::min_element(fits.begin(), fits.end(), [] (auto& a, auto& b) {return a.fval < b.fval;});
    auto best_index = std::distance(fits.begin(), best);

    auto x1 = best->get_parameter("x1").value;
    auto x2 = best->get_parameter("x2").value;

    double w1 = x1;
    double w2 = x2*(1-w1); // w2 = w2*(1-w1)
    double w3 = 1 - w1 - w2; // w3 = (1-w1)*(1-w2)
    std::cout << "Best fit: (w1, w2, w3) = (" << w1 << ", " << w2 << ", " << w3 << ") with value "<< best->fval << std::endl;

    auto[h, k, l] = indices[best_index];
    auto res_his = histograms1[h]*w1 + histograms2[k]*w2 + histograms3[l]*w3;
    auto res_fit = fitter::LinearFitter(dataset, res_his).fit();
    res_fit->add_parameter(mini::FittedParameter("w1", w1, 0));
    res_fit->add_parameter(mini::FittedParameter("w2", w2, 0));
    res_fit->add_parameter(mini::FittedParameter("w3", w3, 0));

    fitter::FitReporter::report(res_fit);
    fitter::FitReporter::save(res_fit, settings::general::output + "fit.txt");
    plots::PlotIntensityFit(res_fit).save(settings::general::output + "fit." + settings::plots::format);
    plots::PlotIntensityFitResiduals::quick_plot(res_fit, settings::general::output + "residuals." + settings::plots::format);

    return 0;
}