#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/Constants.h>
#include <dataset/SimpleDataset.h>

const auto& q_vals = constants::axes::q_vals;
TEST_CASE("ExvFormFactor::evaluate") {}

#include <plots/PlotDataset.h>
TEST_CASE("ExvFormFactor::plot", "[manual]") {
    plots::PlotDataset plot;
    {
        const form_factor::FormFactor& ff_exv = form_factor::storage::get_form_factor(form_factor::form_factor_t::EXCLUDED_VOLUME);
        SimpleDataset dataset;
        dataset.add_plot_options({{"legend", "average exv"}, {"xlabel", "q"}, {"ylabel", "Amplitude"}, {"logx", true}, {"color", style::color::next()}, {"linewidth", 2}});
        for (const double& q : q_vals) {
            dataset.push_back(q, ff_exv.evaluate(q));
        }
        plot.plot(dataset);
    }

    for (unsigned int ff = 0; ff < form_factor::get_count_without_excluded_volume(); ++ff) {
        const form_factor::ExvFormFactor& ff_obj = form_factor::storage::exv::get_form_factor(static_cast<form_factor::form_factor_t>(ff));
        SimpleDataset dataset;
        dataset.add_plot_options({{"legend", form_factor::to_string(static_cast<form_factor::form_factor_t>(ff))}, {"color", style::color::next()}});
        for (const double& q : q_vals) {
            dataset.push_back(q, ff_obj.evaluate_normalized(q));
        }
        plot.plot(dataset);
    }
    plot.save("temp/test/form_factor/exv_form_factors.png");
}