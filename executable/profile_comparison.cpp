#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/detail/XVGReader.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

// calculate the exact aa profile without form factors
[[maybe_unused]] auto exact_aa_debye = [] (const data::Molecule& molecule) {
    container::Container2D<double> distances(molecule.get_atoms().size(), molecule.get_atoms().size());
    auto atoms = molecule.get_atoms();
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = 0; j < atoms.size(); ++j) {
            distances(i, j) = atoms[i].distance(atoms[j]);
        }
    }

    auto ff = form_factor::storage::atomic::C;
    auto qaxis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    auto q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
    hist::ScatteringProfile I(qaxis);
    for (unsigned int q = q0; q < q0+qaxis.bins; ++q) {
        double sum = 0;
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                double qd = constants::axes::q_vals[q]*distances(i, j);
                if (qd < 1e-6) {
                    sum += 1;
                } else {
                    sum += std::sin(qd)/(qd);
                }
            }
        }
        I.index(q-q0) = sum*std::pow(ff.evaluate(constants::axes::q_vals[q]), 2);
    }
    return I;
};

int main(int argc, char const *argv[]) {
    io::ExistingFile pdb;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/saxs_fitter/")->group("General options");
    CLI11_PARSE(app, argc, argv);

    //### GENERATE INTERNAL PLOT ###//
    settings::general::output += pdb.stem() + "/";
    settings::axes::qmin = 1e-4;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = true;
    settings::molecule::use_occupancy = false;
    settings::general::threads = 1;
    settings::grid::min_bins = 20;

    SimpleDataset 
        ausaxs_crysol_aa, ausaxs_crysol_xx, 
        ausaxs_foxs_aa, ausaxs_foxs_xx, 
        ausaxs_pepsi_aa, ausaxs_pepsi_xx,
        ausaxs_aa, ausaxs_xx, ausaxs_ww,
        ausaxs_xx_grid;

    {   // crysol mimic
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::CrysolManager;
        data::Molecule molecule(pdb);
        molecule.add_implicit_hydrogens();

        auto hist = molecule.get_histogram();
        ausaxs_crysol_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_crysol_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }
    {   // foxs mimic
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::FoXSManager;
        data::Molecule molecule(pdb);
        molecule.add_implicit_hydrogens();

        auto hist = molecule.get_histogram();
        ausaxs_foxs_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_foxs_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }
    {   // pepsi mimic
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PepsiManager;
        data::Molecule molecule(pdb);
        molecule.add_implicit_hydrogens();

        auto hist = molecule.get_histogram();
        ausaxs_pepsi_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_pepsi_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }
    {   // ausaxs
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;
        data::Molecule molecule(pdb);
        molecule.add_implicit_hydrogens();
        molecule.generate_new_hydration();

        auto hist = molecule.get_histogram();
        ausaxs_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
        ausaxs_ww = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_ww().as_dataset();
    }
    {
        // ausaxs grid
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface;
        data::Molecule molecule(pdb);
        molecule.add_implicit_hydrogens();
        molecule.generate_new_hydration();

        auto hist = molecule.get_histogram();
        ausaxs_xx_grid = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();        
    }

    //### CHECK FOR PRESENCE OF EXTERNAL DATA IN OUTPUT FOLDER ###//
    io::File crysol(settings::general::output + "crysol.dat");
    io::File foxs(settings::general::output + "foxs.dat");
    io::File pepsi(settings::general::output + "pepsi.dat");
    io::File waxsis(settings::general::output + "waxsresults/waxs_contrib.xvg");

    if (!crysol.exists() || !foxs.exists() || !pepsi.exists()) {
        std::cout << "External data not found in output folder. Skipping external comparison." << std::endl;
        return 0;
    }

    // generate plot
    SimpleDataset 
        crysol_data_xx, crysol_data_aa, crysol_data_ww,
        pepsi_data_xx, pepsi_data_aa, pepsi_data_ww,
        foxs_data_xx, foxs_data_aa, foxs_data_ww,
        waxsis_data_xx, waxsis_data_aa, waxsis_data_ww;
    SimpleDataset waxsis_0, waxsis_1, waxsis_2, waxsis_3, waxsis_4, waxsis_5, waxsis_6, waxsis_7; 
    {
        Dataset tmp(crysol);
        crysol_data_aa = SimpleDataset(tmp.x(), tmp.col(2));
        crysol_data_xx = SimpleDataset(tmp.x(), tmp.col(3));
        crysol_data_ww = SimpleDataset(tmp.x(), tmp.col(4));

        tmp = Dataset(pepsi);
        pepsi_data_aa = SimpleDataset(tmp.x(), tmp.col(2));
        pepsi_data_xx = SimpleDataset(tmp.x(), tmp.col(3));
        pepsi_data_ww = SimpleDataset(tmp.x(), tmp.col(4));

        tmp = Dataset(settings::general::output + "foxs.dat");
        foxs_data_aa = SimpleDataset(tmp.x(), tmp.col(1));
        foxs_data_xx = SimpleDataset(tmp.x(), tmp.col(2));
        foxs_data_ww = SimpleDataset(tmp.x(), tmp.col(3));
    } if (waxsis.exists()) {
        auto tmp = detail::XVGReader::construct_multifile(waxsis);
        waxsis_data_aa = *tmp[0];
        waxsis_data_xx = *tmp[1];

        io::File ww(settings::general::output + "waxsis_ww.dat");
        if (ww.exists()) {
            waxsis_data_ww = SimpleDataset(ww);
        }

        waxsis_data_aa.normalize();
        waxsis_data_xx.normalize();
        waxsis_data_ww.normalize();
    }

    crysol_data_aa.save(settings::general::output + "crysol_aa.dat" );
    foxs_data_aa.save(  settings::general::output + "foxs_aa.dat"   );
    pepsi_data_aa.save( settings::general::output + "pepsi_aa.dat"  );
    ausaxs_aa.save(     settings::general::output + "ausaxs_aa.dat" );

    crysol_data_aa.normalize();
    crysol_data_xx.normalize();
    crysol_data_ww.normalize();
    foxs_data_aa.normalize();
    foxs_data_xx.normalize();
    foxs_data_ww.normalize();
    pepsi_data_aa.normalize();
    pepsi_data_xx.normalize();
    pepsi_data_ww.normalize();
    ausaxs_crysol_aa.normalize();
    ausaxs_crysol_xx.normalize();
    ausaxs_foxs_aa.normalize();
    ausaxs_foxs_xx.normalize();
    ausaxs_pepsi_aa.normalize();
    ausaxs_pepsi_xx.normalize();
    ausaxs_xx.normalize();
    ausaxs_aa.normalize();
    ausaxs_ww.normalize();
    ausaxs_xx_grid.normalize();

    // auto exact_aa = exact_aa_debye(data::Molecule(pdb)).as_dataset();
    // exact_aa.normalize();

    plots::PlotIntensity()
        .plot(crysol_data_aa, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::cyan}, {"title", pdb.stem() + " $I_{aa}$ profiles"}, {"xrange", Limit(1e-2, 1)}}))
        .plot(foxs_data_aa, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}}))
        .plot(pepsi_data_aa, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}}))
        // .plot(waxsis_data_aa, plots::PlotOptions({{"legend", "WAXSiS"}, {"color", style::color::green}}))
        .plot(ausaxs_aa, plots::PlotOptions({{"legend", "AUSAXS"}, {"color", style::color::black}}))
        // .plot(exact_aa, plots::PlotOptions({{"legend", "Exact"}, {"color", style::color::red}, {"linestyle", style::line::dashed}}))
    .save(settings::general::output + "profiles_aa.png");

    plots::PlotIntensity()
        .plot(crysol_data_ww, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::cyan}, {"title", pdb.stem() + " $I_{ww}$ profiles"}, {"xrange", Limit(1e-2, 1)}, {"yrange", Limit(1e-5, 1.1)}}))
        .plot(foxs_data_ww, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}}))
        .plot(pepsi_data_ww, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}}))
        // .plot(waxsis_data_ww, plots::PlotOptions({{"legend", "WAXSiS"}, {"color", style::color::green}}))
        .plot(ausaxs_ww, plots::PlotOptions({{"legend", "AUSAXS"}, {"color", style::color::black}}))
    .save(settings::general::output + "profiles_ww.png");

    plots::PlotIntensity()
        .plot(crysol_data_xx, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::cyan}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"xrange", Limit(1e-2, 1)}}))
        .plot(foxs_data_xx, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}}))
        .plot(pepsi_data_xx, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}}))
        // .plot(waxsis_data_xx, plots::PlotOptions({{"legend", "WAXSiS_xx"}, {"color", style::color::brown}}))
        // .plot(ausaxs_crysol_xx, plots::PlotOptions({{"color", style::color::black}, {"linestyle", style::line::dashed}, {"linewidth", 0.5}}))
        // .plot(ausaxs_foxs_xx, plots::PlotOptions({{"color", style::color::black}, {"linestyle", style::line::dashed}, {"linewidth", 0.5}}))
        // .plot(ausaxs_pepsi_xx, plots::PlotOptions({{"legend", "AUSAXS$_{mimics}$"}, {"color", style::color::black}, {"linestyle", style::line::dashed}, {"linewidth", 0.5}}))
        .plot(ausaxs_xx, plots::PlotOptions({{"legend", "AUSAXS$_{Fraser}$"}, {"color", style::color::black}}))
        .plot(ausaxs_xx_grid, plots::PlotOptions({{"legend", "AUSAXS$_{grid}$"}, {"color", style::color::black}, {"linestyle", style::line::dotted}}))
    .save(settings::general::output + "profiles_xx.png");

    SimpleDataset crysol_diff_xx = crysol_data_xx, crysol_diff_aa = crysol_data_aa;
    for (size_t i = 0; i < crysol_diff_xx.size(); ++i) {
        crysol_diff_xx.y(i) /= ausaxs_crysol_xx.interpolate_x(crysol_diff_xx.x(i), 1);
        crysol_diff_aa.y(i) /= ausaxs_aa.interpolate_x(crysol_diff_aa.x(i), 1);
    }

    SimpleDataset foxs_diff_xx = foxs_data_xx, foxs_diff_aa = foxs_data_aa;
    for (size_t i = 0; i < foxs_diff_xx.size(); ++i) {
        foxs_diff_xx.y(i) /= ausaxs_foxs_xx.interpolate_x(foxs_diff_xx.x(i), 1);
        foxs_diff_aa.y(i) /= ausaxs_foxs_aa.interpolate_x(foxs_diff_aa.x(i), 1);
    }

    SimpleDataset pepsi_diff_xx = pepsi_data_xx, pepsi_diff_aa = pepsi_data_aa, pepsi_diff_ww = pepsi_data_ww;
    for (size_t i = 0; i < pepsi_diff_xx.size(); ++i) {
        pepsi_diff_xx.y(i) /= ausaxs_pepsi_xx.interpolate_x(pepsi_diff_xx.x(i), 1);
        pepsi_diff_aa.y(i) /= ausaxs_pepsi_aa.interpolate_x(pepsi_diff_aa.x(i), 1);
    }

    // SimpleDataset exact_aa_diff = exact_aa;
    // for (size_t i = 0; i < exact_aa_diff.size(); ++i) {
    //     exact_aa_diff.y(i) /= ausaxs_aa.interpolate_x(exact_aa_diff.x(i), 1);
    // }

    plots::PlotDataset()
        .plot(crysol_diff_aa, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "THEIRS / AUSAXS"}, {"title", pdb.stem() + " $I_{aa}$ profiles"}, {"yrange", Limit(0.9, 1.1)}, {"xrange", Limit(1e-2, 1)}, {"color", style::color::cyan}, {"normalize", true}}))
        .plot(foxs_diff_aa, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}, {"normalize", true}}))
        .plot(pepsi_diff_aa, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}, {"normalize", true}}))
        // .plot(exact_aa_diff, plots::PlotOptions({{"legend", "Exact"}, {"color", style::color::red}, {"linestyle", style::line::dashed}, {"normalize", true}}))
        .hline(1, plots::PlotOptions({{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(settings::general::output + "profiles_aa_diff.png");

    plots::PlotDataset()
        .plot(crysol_diff_xx, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "THEIRS / AUSAXS"}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"yrange", Limit(0.9, 1.1)}, {"xrange", Limit(1e-2, 1)}, {"color", style::color::cyan}, {"normalize", true}}))
        .plot(foxs_diff_xx, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}, {"normalize", true}}))
        .plot(pepsi_diff_xx, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}, {"normalize", true}}))
        .hline(1, plots::PlotOptions({{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(settings::general::output + "profiles_xx_diff.png");

    for (size_t i = 0; i < crysol_diff_xx.size(); ++i) {
        crysol_diff_xx.y(i) = crysol_data_xx.y(i) / ausaxs_xx_grid.interpolate_x(crysol_diff_xx.x(i), 1);
    }

    for (size_t i = 0; i < foxs_diff_xx.size(); ++i) {
        foxs_diff_xx.y(i) = foxs_data_xx.y(i) / ausaxs_xx_grid.interpolate_x(foxs_diff_xx.x(i), 1);
    }

    for (size_t i = 0; i < pepsi_diff_xx.size(); ++i) {
        pepsi_diff_xx.y(i) = pepsi_data_xx.y(i) / ausaxs_xx_grid.interpolate_x(pepsi_diff_xx.x(i), 1);
    }

    for (size_t i = 0; i < ausaxs_xx_grid.size(); ++i) {
        ausaxs_xx.y(i) /= ausaxs_xx_grid.y(i);
    }
    
    plots::PlotDataset()
        .plot(crysol_diff_xx, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "Deviation relative to grid"}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"yrange", Limit(0.1, 100)}, {"xrange", Limit(1e-2, 1)}, {"color", style::color::cyan}, {"normalize", true}, {"logy", true}}))
        .plot(foxs_diff_xx, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}, {"normalize", true}}))
        .plot(pepsi_diff_xx, plots::PlotOptions({{"legend", "Pepsi-SAXS"}, {"color", style::color::blue}, {"normalize", true}}))
        .plot(ausaxs_xx, plots::PlotOptions({{"legend", "AUSAXS"}, {"color", style::color::black}, {"normalize", true}}))
        .hline(1, plots::PlotOptions({{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(settings::general::output + "profiles_xx_diff_grid.png");

    crysol_data_xx.save(settings::general::output + "profile_comparison/crysol_xx.dat");
    crysol_data_aa.save(settings::general::output + "profile_comparison/crysol_aa.dat");
    foxs_data_xx.save(settings::general::output + "profile_comparison/foxs_xx.dat");
    foxs_data_aa.save(settings::general::output + "profile_comparison/foxs_aa.dat");
    pepsi_data_xx.save(settings::general::output + "profile_comparison/pepsi_xx.dat");
    pepsi_data_aa.save(settings::general::output + "profile_comparison/pepsi_aa.dat");
    crysol_diff_aa.save(settings::general::output + "profile_comparison/crysol_aa_diff.dat");
    crysol_diff_xx.save(settings::general::output + "profile_comparison/crysol_xx_diff.dat");
    foxs_diff_aa.save(settings::general::output + "profile_comparison/foxs_aa_diff.dat");
    foxs_diff_xx.save(settings::general::output + "profile_comparison/foxs_xx_diff.dat");
    pepsi_diff_aa.save(settings::general::output + "profile_comparison/pepsi_aa_diff.dat");
    pepsi_diff_xx.save(settings::general::output + "profile_comparison/pepsi_xx_diff.dat");
    ausaxs_aa.save(settings::general::output + "profile_comparison/ausaxs_aa.dat");
    ausaxs_xx.save(settings::general::output + "profile_comparison/ausaxs_xx.dat");
    ausaxs_crysol_xx.save(settings::general::output + "profile_comparison/ausaxs_crysol_xx.dat");
    ausaxs_pepsi_xx.save(settings::general::output + "profile_comparison/ausaxs_pepsi_xx.dat");
    ausaxs_foxs_xx.save(settings::general::output + "profile_comparison/ausaxs_foxs_xx.dat");
}