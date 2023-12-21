#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <hydrate/Grid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

#include <sstream>

int main(int argc, char const *argv[]) {
    std::string s_pdb;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/exv_comparison/")->group("General options");
    CLI11_PARSE(app, argc, argv);

    //### GENERATE INTERNAL PLOT ###//
    io::ExistingFile pdb(s_pdb);
    settings::general::output += pdb.stem() + "/";
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;

    settings::grid::exv_radius = 0.5;
    hist::CompositeDistanceHistogramFFGrid::regenerate_table();
    data::Molecule molecule(pdb);

    auto mtffg1  = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFGrid    <true>(molecule).calculate_all().get())->get_profile_xx().as_dataset();
    auto mtffavg = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFAvg     <true>(molecule).calculate_all().get())->get_profile_xx().as_dataset();
    auto mtffexp = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFExplicit<true>(molecule).calculate_all().get())->get_profile_xx().as_dataset();

    settings::grid::exv_radius = 1;
    molecule.clear_grid();
    hist::CompositeDistanceHistogramFFGrid::regenerate_table();
    auto mtffg2   = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFGrid<true>(molecule).calculate_all().get())->get_profile_xx().as_dataset();

    settings::grid::exv_radius = 1.5;
    molecule.clear_grid();
    hist::CompositeDistanceHistogramFFGrid::regenerate_table();
    auto mtffg3   = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFGrid<true>(molecule).calculate_all().get())->get_profile_xx().as_dataset();

    double cg1 = mtffg1.normalize();
    double cg2 = mtffg2.normalize();
    double cg3 = mtffg3.normalize();
    double ca  = mtffavg.normalize();
    double ce  = mtffexp.normalize();

    std::stringstream ssg1, ssa, sse, ssg2, ssg3;
    ssg1 << "Grid-based (1), c=" << std::fixed << std::setprecision(2) << 1;
    ssg2 << "Grid-based (2), c=" << std::fixed << std::setprecision(2) << cg2/cg1;
    ssg3 << "Grid-based (3), c=" << std::fixed << std::setprecision(2) << cg3/cg1;
    ssa << "Average, c=" << std::fixed << std::setprecision(2) << ca/cg1;
    sse << "Explicit, c=" << std::fixed << std::setprecision(2) << ce/cg1;

    plots::PlotIntensity()
        .plot(mtffg1, plots::PlotOptions({{"legend", ssg1.str()}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::next()}, {"title", pdb.stem()}}))
        .plot(mtffg2, plots::PlotOptions({{"legend", ssg2.str()}, {"color", style::color::next()}}))
        .plot(mtffg3, plots::PlotOptions({{"legend", ssg3.str()}, {"color", style::color::next()}}))
        .plot(mtffavg, plots::PlotOptions({{"legend", ssa.str()}, {"color", style::color::next()}}))
        .plot(mtffexp, plots::PlotOptions({{"legend", sse.str()}, {"color", style::color::next()}}))
    .save(settings::general::output + "intensity.png");


    //### CHECK FOR PRESENCE OF EXTERNAL DATA IN OUTPUT FOLDER ###//
    io::File crysol(settings::general::output + "crysol.dat");
    io::File foxs(settings::general::output + "foxs_xx.dat");
    io::File gromacs_prot_pdb(settings::general::output + "prot+solvlayer_0.pdb");
    io::File gromacs_exv_pdb(settings::general::output + "excludedvolume_0.pdb");

    if (!crysol.exists() || !foxs.exists() || !gromacs_prot_pdb.exists() || !gromacs_exv_pdb.exists()) {
        std::cout << "External data not found in output folder. Skipping external comparison." << std::endl;
        return 0;
    }

    // generate plot
    SimpleDataset crysol_data;
    {
        Dataset tmp(crysol);
        crysol_data = SimpleDataset(tmp.x(), tmp.col(3));
    }
    auto foxs_data = SimpleDataset(foxs);

    settings::molecule::use_effective_charge = false;
    settings::molecule::throw_on_unknown_atom = false;
    SimpleDataset gromacs_data;
    {
        data::Molecule gromacs_x(gromacs_exv_pdb);
        data::Molecule gromacs_p(gromacs_prot_pdb);
        gromacs_p.clear_hydration();

        std::vector<data::record::Water> waters;
        waters.reserve(gromacs_x.get_waters().size());
        auto grid = gromacs_p.get_grid();
        grid->expand_volume();
        for (auto& w : gromacs_x.get_waters()) {
            auto p = grid->to_bins_bounded(w.get_coordinates());
            if (grid->grid.is_atom_area_or_volume(p.x(), p.y(), p.z())) {
                waters.push_back(w);
            }
        }
        data::Molecule gromacs(std::vector<data::record::Atom>{}, waters);
        gromacs_data = gromacs.get_histogram()->debye_transform().as_dataset();
    }

    crysol_data.normalize();
    foxs_data.normalize();
    gromacs_data.normalize();

    plots::PlotIntensity()
        .plot(crysol_data, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::cyan}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"yrange", Limit(1e-4, 1.1)}, {"xrange", Limit(1e-2, 1)}}))
        .plot(foxs_data, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}}))
        .plot(gromacs_data, plots::PlotOptions({{"legend", "GROMACS"}, {"color", style::color::green}}))
        .plot(mtffg2, plots::PlotOptions({{"legend", "Grid-based"}, {"color", style::color::red}}))
        .plot(mtffavg, plots::PlotOptions({{"legend", "Average"}, {"color", style::color::purple}}))
        .plot(mtffexp, plots::PlotOptions({{"legend", "Explicit"}, {"color", style::color::brown}}))
    .save(settings::general::output + "comparison.png");
}