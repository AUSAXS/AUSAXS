#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/Constants.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>

const auto& q_vals = constants::axes::q_vals;
TEST_CASE("ExvFormFactor::evaluate") {}

// compare the excluded volume form factors with the average one derived by Jan
TEST_CASE("ExvFormFactor::plot", "[manual]") {
    plots::PlotDataset plot;
    {
        const form_factor::FormFactor& ff_exv = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::EXCLUDED_VOLUME);
        SimpleDataset dataset;
        for (const double& q : q_vals) {
            dataset.push_back(q, ff_exv.evaluate(q)/ff_exv.evaluate(0));
        }
        plot.plot(dataset, plots::PlotOptions({{"legend", "average exv"}, {"xlabel", "q"}, {"ylabel", "Amplitude"}, {"logx", true}, {"logy", true}, {"color", style::color::next()}, {"linewidth", 2}}));
    }

    for (unsigned int ff = 0; ff < form_factor::get_count_without_excluded_volume(); ++ff) {
        const form_factor::ExvFormFactor& ff_obj = form_factor::storage::exv::standard.get_form_factor(static_cast<form_factor::form_factor_t>(ff));
        SimpleDataset dataset;
        for (const double& q : q_vals) {
            dataset.push_back(q, ff_obj.evaluate_normalized(q));
        }
        plot.plot(dataset, plots::PlotOptions({{"legend", form_factor::to_string(static_cast<form_factor::form_factor_t>(ff))}, {"color", style::color::next()}}));
    }
    plot.save("temp/tests/form_factor/exv_form_factors.png");
}

// compare each exv form factor with its real one
TEST_CASE("ExvFormFactor::plot_cmp", "[manual]") {
    for (unsigned int ffi = 0; ffi < form_factor::get_count_without_excluded_volume(); ++ffi) {
        const form_factor::FormFactor& ff = form_factor::storage::atomic::get_form_factor(static_cast<form_factor::form_factor_t>(ffi));
        const form_factor::ExvFormFactor& ffx = form_factor::storage::exv::standard.get_form_factor(static_cast<form_factor::form_factor_t>(ffi));

        SimpleDataset dataset, datasetx;
        for (const double& q : q_vals) {
            dataset.push_back(q, ff.evaluate(q)*ffx.evaluate(0));
            datasetx.push_back(q, ffx.evaluate(q));
        }
        plots::PlotDataset()
            .plot(dataset, plots::PlotOptions({{"legend", form_factor::to_string(static_cast<form_factor::form_factor_t>(ffi))}, {"color", style::color::orange}}))
            .plot(datasetx, plots::PlotOptions({{"legend", form_factor::to_string(static_cast<form_factor::form_factor_t>(ffi)) + "x"}, {"color", style::color::black}}))
        .save("temp/tests/form_factor/cmp/" + form_factor::to_string(static_cast<form_factor::form_factor_t>(ffi)) + ".png");
    }
}