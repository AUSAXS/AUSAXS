// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_debye.h>
#include <api/cli/detail/CLICommon.h>

#include <CLI/CLI.hpp>

#include <constants/ValidFileExtensions.h>
#include <data/Molecule.h>
#include <dataset/Dataset.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <hist/scattering_model/ScatteringModelPoint.h>
#include <io/ExistingFile.h>
#include <settings/All.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <string>
#include <vector>

using namespace ausaxs;

int cli_debye(int argc, char const *argv[]) {
    io::ExistingFile structure;
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
    cli::add_general_options(app, "output/debye/");
    cli::add_structure_options(app);
    app.add_option("--qmin", qmin, "Lower limit of the q-range in Å^-1.")->default_val(qmin)->check(CLI::NonNegativeNumber);
    app.add_option("--qmax", qmax, "Upper limit of the q-range in Å^-1.")->default_val(qmax)->check(CLI::PositiveNumber);
    app.add_option("--count,-n", count, "The number of q values, evenly spaced from qmin to qmax.")->default_val(count)->check(CLI::Range(2, 1000000));
    app.add_flag("--points", points, "Treat every atom as a point scatterer with its weight, instead of using its form factor.");
    app.add_option_function<double>("--bin-width", [] (double width) {settings::axes::bin_width = width;}, // through the setting, which also updates the inverse width
        "The width of the distance bins in Å. Lower it for large crystals.")->check(CLI::PositiveNumber);
    auto* p_wavelength = app.add_option("--wavelength", wavelength, "Also write the intensity against the scattering angle 2θ for this wavelength in Å.")->check(CLI::PositiveNumber);
    app.add_flag("--neutron", neutron, "Use neutron scattering lengths instead of X-ray form factors.");
    cli::add_config_subcommand(app, config, "debye");
    cli::add_molecule_subcommand(app);
    cli::add_hydrogens_subcommand(app);
    cli::add_neutron_solvent_options(app.add_subcommand("solv", "The solvent the exchangeable hydrogens take their isotope composition from. Only used with --neutron."));

    app.final_callback([&] () {
        if (config.save) {
            settings::write("settings.txt");
            console::print_info("Settings saved to settings.txt in current directory.");
            if (!input_s->count()) {throw CLI::Success();}
        }

        // required args (not marked ->required() since that interferes with the help flag for subcommands)
        if (!input_s->count()) {
            console::print_warning("Error: input_structure is required.");
            throw CLI::RuntimeError(1);
        }
    });

    try {
        if (auto code = cli::parse(app, argc, argv, config, [&] () {return structure.directory();})) {return *code;}

        // the defaults of this tool, which take priority over a settings file
        settings::scattering::radiation = neutron ? settings::scattering::Radiation::Neutron : settings::scattering::Radiation::XRay;
        settings::exv::exv_method = points ? settings::exv::ExvMethod::None : settings::exv::ExvMethod::Vacuum;
        settings::hist::weighted_bins = true; // the profile must hold far beyond the range unweighted bins are judged by
        settings::validate_settings();

        if (!constants::filetypes::structure.check(structure)) {
            throw except::invalid_argument("Unknown structure extension: " + structure.str());
        }
        if (qmax <= qmin) {
            throw except::invalid_argument("qmax (" + std::to_string(qmax) + ") must be larger than qmin (" + std::to_string(qmin) + ").");
        }

        if (neutron) {cli::print_neutron_solvent();}
        data::Molecule molecule(structure);
        if (molecule.size_water() != 0) {
            console::print_text("\tDiscarding " + std::to_string(molecule.size_water()) + " water molecules from the structure file.");
            molecule.clear_hydration();
        }

        std::vector<double> q(count);
        for (int i = 0; i < count; ++i) {q[i] = qmin + (qmax - qmin)*i/(count - 1);}
        Dataset I = points
            ? hist::ScatteringModelPoint(std::move(*molecule.get_total_histogram())).debye_transform(q)
            : molecule.get_histogram()->debye_transform(q)
        ;

        settings::general::output += structure.stem() + "/";
        const std::string description = std::string(points ? "point scatterers" : (neutron ? "neutron scattering lengths" : "X-ray form factors"));
        I.save(settings::general::output + "debye.dat", "q I, " + std::to_string(molecule.size_atom()) + " atoms with " + description);
        console::print_text("Wrote I(q) of " + std::to_string(molecule.size_atom()) + " atoms with " + description + " to " + settings::general::output + "debye.dat");

        // q = 4π/λ sin(θ), so only q ≤ 4π/λ is reachable at this wavelength
        if (p_wavelength->count() != 0) {
            std::vector<double> two_theta, I_theta;
            for (int i = 0; i < I.size(); ++i) {
                double s = I.x(i)*wavelength/(4*std::numbers::pi);
                if (1 < s) {break;}
                two_theta.push_back(2*std::asin(s)*180/std::numbers::pi);
                I_theta.push_back(I.y(i));
            }
            if (static_cast<int>(two_theta.size()) < I.size()) {
                console::print_text("\tq > " + std::to_string(4*std::numbers::pi/wavelength) + " Å^-1 cannot be reached at a wavelength of " + std::to_string(wavelength) + " Å and is left out of the 2θ pattern.");
            }
            // a plain table, since a Dataset annotates its first column as q
            std::ofstream out(settings::general::output + "debye_2theta.dat");
            out << "2theta[deg] I, wavelength " + std::to_string(wavelength) + " A\n\n" << std::scientific << std::setprecision(8);
            for (std::size_t i = 0; i < two_theta.size(); ++i) {out << two_theta[i] << "\t" << I_theta[i] << "\n";}
            console::print_text("Wrote I(2θ) to " + settings::general::output + "debye_2theta.dat");
        }
    } catch (const except::base&) {
        return 1; // our own exceptions already report themselves when thrown
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        return 1;
    }
    return 0;
}
