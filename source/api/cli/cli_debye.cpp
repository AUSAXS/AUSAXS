// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_debye.h>
#include <api/cli/detail/CLICommon.h>

#include <CLI/CLI.hpp>

#include <constants/ValidFileExtensions.h>
#include <data/Molecule.h>
#include <dataset/Dataset.h>
#include <dataset/SimpleDataset.h>
#include <fitter/FitReporter.h>
#include <fitter/FitResult.h>
#include <fitter/detail/LinearLeastSquares.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <hist/scattering_model/ScatteringModelPoint.h>
#include <io/ExistingFile.h>
#include <settings/All.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

using namespace ausaxs;

namespace {
    /**
     * @brief Write columns against the scattering angle 2θ in degrees, for q = 4π/λ sin(θ). 
     *        Only q ≤ 4π/λ can be reached at this wavelength, so larger q are left out.
     *        A plain table, since a Dataset annotates its first column as q.
     */
    void write_two_theta(const std::string& path, double wavelength, const std::vector<double>& q, const std::vector<std::vector<double>>& cols, const std::string& names) {
        std::size_t n = 0;
        while (n < q.size() && q[n]*wavelength/(4*std::numbers::pi) <= 1) {++n;}
        if (n < q.size()) {
            console::print_text("\tq > " + std::to_string(4*std::numbers::pi/wavelength) + " Å^-1 cannot be reached at a wavelength of " + std::to_string(wavelength) + " Å and is left out of the 2θ pattern.");
        }
        std::ofstream out(path);
        out << "2theta[deg] " << names << ", wavelength " << std::to_string(wavelength) << " A\n\n" << std::scientific << std::setprecision(8);
        for (std::size_t i = 0; i < n; ++i) {
            out << 2*std::asin(q[i]*wavelength/(4*std::numbers::pi))*180/std::numbers::pi;
            for (const auto& col : cols) {out << "\t" << col[i];}
            out << "\n";
        }
        console::print_text("Wrote the pattern against 2θ to " + path);
    }
}

int cli_debye(int argc, char const *argv[]) {
    io::ExistingFile structure, mfile;
    cli::ConfigOptions config;
    double qmin = 1e-4, qmax = 1, wavelength = 0;
    int count = 1000;
    bool points = false, neutron = false;

    CLI::App app{
        "Calculate the Debye scattering intensity of a structure, such as a crystal, a nanoparticle or a biomolecule, "
        "without the solvent: there is no excluded volume and no hydration shell."
    };
    app.fallthrough();
    auto* input_s = app.add_option("input_structure", structure, "Path to the structure file.")->check(CLI::ExistingFile);
    auto* input_m = app.add_option("input_measurement", mfile, "Path to a measurement to compare with. The intensity is then evaluated at its q values.")->check(CLI::ExistingFile);
    cli::add_general_options(app, "output/debye/");
    cli::add_structure_options(app);
    // not the q-range settings, which are bound to the predefined q-axis, while this range may extend far beyond it
    auto* p_qmin = app.add_option("--qmin", qmin, "Lower limit of the q-range in Å^-1. With a measurement, it only narrows the measured range.")->default_val(qmin)->check(CLI::NonNegativeNumber);
    auto* p_qmax = app.add_option("--qmax", qmax, "Upper limit of the q-range in Å^-1. With a measurement, it only narrows the measured range.")->default_val(qmax)->check(CLI::PositiveNumber);
    auto* p_count = app.add_option("--count,-n", count, "The number of q values, evenly spaced from qmin to qmax. Not used with a measurement.")->default_val(count)->check(CLI::Range(2, 1000000));
    app.add_flag("--points", points, "Treat every atom as a point scatterer with its weight, instead of using its form factor.");
    app.add_option_function<double>("--bin-width", [] (double width) {settings::axes::bin_width = width;}, // through the setting, which also updates the inverse width
        "The width of the distance bins in Å. Lower it for large crystals.")->check(CLI::PositiveNumber);
    auto* p_wavelength = app.add_option("--wavelength", wavelength, "Also write the intensity against the scattering angle 2θ for this wavelength in Å.")->check(CLI::PositiveNumber);
    app.add_flag("--neutron", neutron, "Use neutron scattering lengths instead of X-ray form factors.");
    cli::add_config_subcommand(app, config, "debye");
    cli::add_molecule_subcommand(app);
    cli::add_hydrogens_subcommand(app);
    cli::add_measurement_options(app.add_subcommand("data", "See and set additional options for reading the measurement."));
    cli::add_neutron_solvent_options(app.add_subcommand("solv", "The solvent the exchangeable hydrogens take their isotope composition from. Only used with --neutron."));

    app.final_callback([&] () {
        if (config.save) {
            settings::write("settings.txt");
            console::print_info("Settings saved to settings.txt in current directory.");
            if (input_s->count() == 0) {throw CLI::Success();}
        }

        // required args (not marked ->required() since that interferes with the help flag for subcommands)
        if (input_s->count() == 0) {
            console::print_warning("Error: input_structure is required.");
            throw CLI::RuntimeError(1);
        }
    });

    try {
        if (auto code = cli::parse(app, argc, argv, config, [&] () {return input_m->count() != 0 ? mfile.directory() : structure.directory();})) {return *code;}

        // the defaults of this tool, which take priority over a settings file
        settings::scattering::radiation = neutron ? settings::scattering::Radiation::Neutron : settings::scattering::Radiation::XRay;
        settings::exv::exv_method = points ? settings::exv::ExvMethod::None : settings::exv::ExvMethod::Vacuum;
        settings::hist::weighted_bins = true; // the profile must hold far beyond the range unweighted bins are judged by
        settings::validate_settings();

        if (!constants::filetypes::structure.check(structure)) {
            throw except::invalid_argument("Unknown structure extension: " + structure.str());
        }
        if ((input_m->count() == 0 || (p_qmin->count() != 0 && p_qmax->count() != 0)) && qmax <= qmin) { // with a measurement, an unset limit is its own
            throw except::invalid_argument("qmax (" + std::to_string(qmax) + ") must be larger than qmin (" + std::to_string(qmin) + ").");
        }
        if (input_m->count() != 0 && !constants::filetypes::saxs_data.check(mfile)) {
            throw except::invalid_argument("Unknown measurement extension: " + mfile.str());
        }

        if (neutron) {cli::print_neutron_solvent();}
        data::Molecule molecule(structure);
        if (molecule.size_water() != 0) {
            console::print_text("\tDiscarding " + std::to_string(molecule.size_water()) + " water molecules from the structure file.");
            molecule.clear_hydration();
        }

        auto debye = [&] (const std::vector<double>& q) -> std::vector<double> {
            auto I = points 
                ? hist::ScatteringModelPoint(std::move(*molecule.get_total_histogram())).debye_transform(q)
                : molecule.get_histogram()->debye_transform(q)
            ;
            return I.y();
        };
        const std::string description = std::string(points ? "point scatterers" : (neutron ? "neutron scattering lengths" : "X-ray form factors"));

        // fitting mode: the scale and background are always fitted
        if (input_m->count() != 0) {
            // the measurement may extend far beyond the q-range settings, so it is read whole and narrowed here instead
            settings::axes::clamp_to_qrange = false;
            settings::internal_state::suggest_rebinning = false; // rebinning would smear the peaks of a diffraction pattern
            SimpleDataset measurement(mfile);
            measurement.limit(
                0, 
                p_qmin->count() != 0 ? qmin : -std::numeric_limits<double>::infinity(), 
                p_qmax->count() != 0 ? qmax : std::numeric_limits<double>::infinity()
            );
            if (measurement.size() < 3) {throw except::invalid_argument("The measurement \"" + mfile.str() + "\" has fewer than three q values in the requested range.");}
            if (p_count->count() != 0) {console::print_text("\tIgnoring --count, since the q values are taken from the measurement.");}
            settings::general::output += mfile.stem() + "/";

            // the model is evaluated directly at the measured q values, which may lie outside the predefined q-axis the SmartFitter works on
            fitter::detail::LinearLeastSquares linear_fitter(debye(measurement.x()), measurement.y(), measurement.yerr());
            auto result = linear_fitter.fit();
            auto model = linear_fitter.get_model_curve(result->get_parameter_values());
            result->set_data_curves(measurement.x(), measurement.y(), measurement.yerr(), std::vector<double>(model), linear_fitter.get_residuals(result->get_parameter_values()));
            fitter::FitReporter::report(result.get());
            fitter::FitReporter::save(result.get(), settings::general::output + "report.txt", argc, argv);
            result->curves.select_columns({0, 1, 2, 3}).save(
                settings::general::output + "ausaxs.fit", 
                "chi2=" + std::to_string(result->fval/result->dof) + " dof=" + std::to_string(result->dof) + ", " + std::to_string(molecule.size_atom()) + " atoms with " + description
            );
            console::print_text("Wrote the fit of " + std::to_string(molecule.size_atom()) + " atoms with " + description + " to " + settings::general::output + "ausaxs.fit");
            if (p_wavelength->count() != 0) {
                write_two_theta(settings::general::output + "ausaxs_2theta.fit", wavelength, measurement.x(), {measurement.y(), measurement.yerr(), model}, "I I_err I_fit");
            }
        }

        // simulation mode: the raw profile
        else {
            console::print_info("\nSimulation mode enabled.");
            std::vector<double> q(count);
            for (int i = 0; i < count; ++i) {q[i] = qmin + (qmax - qmin)*i/(count - 1);}
            auto I = debye(q);

            settings::general::output += "simulated/" + structure.stem() + "/";
            Dataset({q, I}).save(settings::general::output + "debye.dat", "q I, " + std::to_string(molecule.size_atom()) + " atoms with " + description);
            console::print_text("Wrote I(q) of " + std::to_string(molecule.size_atom()) + " atoms with " + description + " to " + settings::general::output + "debye.dat");
            if (p_wavelength->count() != 0) {
                write_two_theta(settings::general::output + "debye_2theta.dat", wavelength, q, {I}, "I");
            }
        }
    } catch (const except::base&) {
        return 1; // our own exceptions already report themselves when thrown
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        return 1;
    }
    return 0;
}
