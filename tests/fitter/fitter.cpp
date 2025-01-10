#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <em/ImageStack.h>
#include <plots/All.h>
#include <fitter/FitReporter.h>
#include <utility/Utility.h>
#include <settings/All.h>
#include <hist/Histogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <em/detail/ExtendedLandscape.h>
#include <fitter/SmartFitter.h>
#include <mini/detail/FittedParameter.h>

using namespace ausaxs;
using namespace data;

// TEST_CASE("fitter: consistent_charge_scaling") {
//     settings::molecule::use_effective_charge = true;
//     settings::general::verbose = false;
//     settings::hist::histogram_manager = GENERATE(
//         settings::hist::HistogramManagerChoice::HistogramManager, 
//         settings::hist::HistogramManagerChoice::PartialHistogramManager,
//         settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
//         settings::hist::HistogramManagerChoice::HistogramManagerMT,
//         settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg,
//         settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg,
//         settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit,
//         settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid
//     );

//     auto string_mapper = [] (settings::hist::HistogramManagerChoice choice) -> std::string {
//         switch (choice) {
//             case settings::hist::HistogramManagerChoice::HistogramManager: return "HistogramManager";
//             case settings::hist::HistogramManagerChoice::PartialHistogramManager: return "PartialHistogramManager";
//             case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT: return "PartialHistogramManagerMT";
//             case settings::hist::HistogramManagerChoice::HistogramManagerMT: return "HistogramManagerMT";
//             case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg: return "HistogramManagerMTFFAvg";
//             case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit: return "HistogramManagerMTFFExplicit";
//             case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: return "HistogramManagerMTFFGrid";
//             default: return "Unknown";
//         }
//     };

//     SECTION("simple, " + string_mapper(settings::hist::histogram_manager)) {
//         io::ExistingFile mfile = "tests/files/2epe.dat";
//         Molecule protein("tests/files/2epe.pdb");
//         protein.clear_hydration();
//         fitter::SmartFitter fitter(mfile, protein.get_histogram());
//         auto fit1 = fitter.fit_chi2_only();

//         protein.update_effective_charge(1.2);
//         fitter.set_model(protein.get_histogram());
//         auto fit2 = fitter.fit_chi2_only();
//         REQUIRE(fit1 != fit2);

//         protein.update_effective_charge(1.0);
//         fitter.set_model(protein.get_histogram());
//         fit2 = fitter.fit_chi2_only();
//         REQUIRE_THAT(fit1, Catch::Matchers::WithinAbs(fit2, 1e-6));
//     }
// }

class SmartFitterDebug : public fitter::SmartFitter {
    public:
        using fitter::SmartFitter::SmartFitter;
        void set_data(SimpleDataset& data) {this->data = data;}
};

TEST_CASE("SmartFitter::fit") {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;
    Molecule protein("tests/files/2epe.pdb");

    SmartFitterDebug fitter({{}}, protein.get_histogram());
    auto h = static_cast<hist::ICompositeDistanceHistogramExv*>(fitter.get_model());

    SECTION("hydration shell") {
        settings::fit::fit_hydration = true;
        settings::fit::fit_excluded_volume = false;
        settings::fit::fit_solvent_density = false;

        for (auto v : std::vector{0.5, 2., 5.}) {
            h->apply_water_scaling_factor(v);
            auto d = h->debye_transform().as_dataset();
            d.simulate_errors();
            fitter.set_data(d);
            REQUIRE_THAT(fitter.fit()->get_parameter(constants::fit::Parameters::SCALING_WATER).value, Catch::Matchers::WithinAbs(v, 1e-3));
        }
    }

    SECTION("excluded volume") {
        settings::fit::fit_hydration = false;
        settings::fit::fit_excluded_volume = true;
        settings::fit::fit_solvent_density = false;

        for (auto v : std::vector{0.92, 0.96, 1.04, 1.08}) {
            h->apply_excluded_volume_scaling_factor(v);
            auto d = h->debye_transform().as_dataset();
            d.simulate_errors();
            fitter.set_data(d);
            REQUIRE_THAT(fitter.fit()->get_parameter(constants::fit::Parameters::SCALING_EXV).value, Catch::Matchers::WithinAbs(v, 1e-3));
        }
    }

    SECTION("solvent density") {
        settings::fit::fit_hydration = false;
        settings::fit::fit_excluded_volume = false;
        settings::fit::fit_solvent_density = true;

        for (auto v : std::vector{0.95, 0.975, 1.025, 1.05}) {
            h->apply_solvent_density_scaling_factor(v);
            auto d = h->debye_transform().as_dataset();
            d.simulate_errors();
            fitter.set_data(d);
            REQUIRE_THAT(fitter.fit()->get_parameter(constants::fit::Parameters::SCALING_RHO).value, Catch::Matchers::WithinAbs(v, 1e-3));
        }
    }

    SECTION("atomic debye waller") {
        settings::fit::fit_hydration = false;
        settings::fit::fit_excluded_volume = false;
        settings::fit::fit_solvent_density = false;
        settings::fit::fit_atomic_debye_waller = true;

        for (auto v : std::vector{0.1, 0.5, 1.0, 2.0}) {
            h->apply_atomic_debye_waller_factor(v);
            auto d = h->debye_transform().as_dataset();
            d.simulate_errors();
            fitter.set_data(d);
            REQUIRE_THAT(fitter.fit()->get_parameter(constants::fit::Parameters::DEBYE_WALLER_ATOMIC).value, Catch::Matchers::WithinAbs(v, 1e-3));
        }
    }

    SECTION("excluded volume debye waller") {
        settings::fit::fit_hydration = false;
        settings::fit::fit_excluded_volume = false;
        settings::fit::fit_solvent_density = false;
        settings::fit::fit_atomic_debye_waller = false;
        settings::fit::fit_exv_debye_waller = true;

        for (auto v : std::vector{0.1, 0.5, 1.0, 2.0}) {
            h->apply_exv_debye_waller_factor(v);
            auto d = h->debye_transform().as_dataset();
            d.simulate_errors();
            fitter.set_data(d);
            REQUIRE_THAT(fitter.fit()->get_parameter(constants::fit::Parameters::DEBYE_WALLER_EXV).value, Catch::Matchers::WithinAbs(v, 1e-3));
        }
    }
}

TEST_CASE("fitter: correct dof") {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;
    Molecule protein("tests/files/2epe.pdb");
    SimpleDataset data("tests/files/2epe.dat");
    unsigned int size = data.size();

    settings::fit::fit_hydration = false;
    settings::fit::fit_excluded_volume = false;
    settings::fit::fit_solvent_density = false;
    settings::fit::fit_atomic_debye_waller = false;
    settings::fit::fit_exv_debye_waller = false;

    SECTION("LinearFitter") {
        fitter::LinearFitter fitter(data, protein.get_histogram());
        REQUIRE(fitter.dof() == size-2);
        auto res = fitter.fit();
        REQUIRE(res->dof == size-2);
    }

    SECTION("SmartFitter") {
        auto fit_and_check = [size, &protein, &data] (unsigned int dof) {
            fitter::SmartFitter fitter(data);
            REQUIRE(fitter.dof() == size-dof);
            fitter.set_model(protein.get_histogram());
            auto res = fitter.fit();
            REQUIRE(res->dof == size-dof);
        };

        SECTION("single") {
            settings::fit::fit_hydration = true;
            fit_and_check(3);

            settings::fit::fit_hydration = false;
            settings::fit::fit_excluded_volume = true;
            fit_and_check(3);

            settings::fit::fit_excluded_volume = false;
            settings::fit::fit_solvent_density = true;
            fit_and_check(3);

            settings::fit::fit_solvent_density = false;
            settings::fit::fit_atomic_debye_waller = true;
            fit_and_check(3);

            settings::fit::fit_atomic_debye_waller = false;
            settings::fit::fit_exv_debye_waller = true;
            fit_and_check(3);
        }

        SECTION("double") {
            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = false;
            fit_and_check(4);

            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = false;
            settings::fit::fit_solvent_density = true;
            fit_and_check(4);

            settings::fit::fit_hydration = false;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = true;
            fit_and_check(4);

            settings::fit::fit_hydration = false;
            settings::fit::fit_excluded_volume = false;
            settings::fit::fit_solvent_density = false;
            settings::fit::fit_atomic_debye_waller = true;
            settings::fit::fit_exv_debye_waller = true;
            fit_and_check(4);
        }

        SECTION("triple") {
            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = true;
            fit_and_check(5);
        }

        SECTION("quadruple") {
            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = true;
            settings::fit::fit_atomic_debye_waller = true;
            fit_and_check(6);

            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = true;
            settings::fit::fit_atomic_debye_waller = false;
            settings::fit::fit_exv_debye_waller = true;
            fit_and_check(6);
        }

        SECTION("quintuple") {
            settings::fit::fit_hydration = true;
            settings::fit::fit_excluded_volume = true;
            settings::fit::fit_solvent_density = true;
            settings::fit::fit_atomic_debye_waller = true;
            settings::fit::fit_exv_debye_waller = true;
            fit_and_check(7);
        }
    }
}