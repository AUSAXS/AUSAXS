#include <CLI/CLI.hpp>

#include <grid/Grid.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <constants/Constants.h>
#include <table/ArrayDebyeTable.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hydrate/generation/RadialHydration.h>
#include <dataset/SimpleDataset.h>
#include <form_factor/FormFactor.h>
#include <settings/All.h>
#include <plots/All.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);

    io::Folder input;
    io::File file_unordered, file_ordered;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("folder", input, "Path to the MD SAXS folder.")->required()->check(CLI::ExistingDirectory);
    app.add_option("--u", file_unordered, "Path to the unordered envelope file.")->default_val("excludedvolume_0.pdb");
    app.add_option("--o", file_ordered, "Path to the ordered envelope file.")->default_val("prot+solvlayer_0.pdb");
    CLI11_PARSE(app, argc, argv);

    std::string name = input.path();
    name = name.substr(name.find_last_of("_")+1);
    {
        auto nested = input.directories();
        if (nested.size() != 1) {
            std::cerr << "Error: The input folder must contain exactly one subfolder." << std::endl;
            return 1;
        }
        input = nested[0];
        nested = input.directories();
        if (nested.size() != 1) {
            std::cerr << "Error: The nested input folder must contain exactly one subfolder." << std::endl;
            return 1;
        }
        if (!nested.front().path().ends_with("envelope")) {
            std::cout << nested.front().path() << std::endl;
            std::cerr << "Error: The input folder must contain a subfolder called 'envelope'." << std::endl;
            return 1;
        }
        input = nested.front();
    }

    bool calc_scattering = false;
    bool calc_density = true;

    // settings::hydrate::shell_correction = 0;
    // hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0,0,0};});

    settings::grid::scaling = 1;
    settings::grid::cell_width = 0.5;
    settings::grid::min_exv_radius = 2.5; // water + atom
    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::general::output = "output/waxsis_water/";
    settings::general::keep_hydrogens = false;

    // io::ExistingFile env_unordered(input + file_unordered.filename());
    io::ExistingFile env_ordered(input + file_ordered.filename());

    // env_unordered = env_unordered.copy("temp/ausaxs");
    env_ordered = env_ordered.copy("temp/ausaxs");

    auto remove_cip = [] (std::string path) {
        std::ifstream in(path);
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(in, line)) {
            if (line.find("CIP") == std::string::npos && line.find("CIM") == std::string::npos) {
                lines.push_back(line);
            }
        }

        std::ofstream out(path);
        for (auto& line : lines) {
            out << line << std::endl;
        }
    };
    remove_cip(env_ordered);

    data::Molecule ordered(env_ordered);
    // data::Molecule disordered(env_unordered);
    // disordered.get_body(0).get_atoms() = ordered.get_body(0).get_atoms();
    ordered.save(settings::general::output + "ordered.pdb");
    // disordered.save(settings::general::output + "disordered.pdb");

    auto impose_exv = [] (data::Molecule& protein) {
        auto old_waters = protein.get_waters();
        protein.clear_hydration();
        auto grid = protein.get_grid();
        grid->expand_volume();

        std::vector<data::record::Water> new_waters;
        for (auto& water : old_waters) {
            auto bin = grid->to_bins(water.get_coordinates());
            if (grid->grid.is_empty(bin.x(), bin.y(), bin.z())) {
                new_waters.emplace_back(water);
            }
        }
        return data::Molecule(protein.get_atoms(), new_waters);
    };

    // ordered = impose_exv(ordered);
    // disordered = impose_exv(disordered);

    //#######################//
    //### CALC SCATTERING ###//
    //#######################//
    if (calc_scattering) {
        auto sinqd_table = table::ArrayDebyeTable::get_default_table();
        std::vector<data::record::Water> waters = ordered.get_waters();
        // waters.insert(waters.end(), disordered.get_waters().begin(), disordered.get_waters().end());
        // for (unsigned int i = disordered.size_water(); i < waters.size(); ++i) {
        //     waters[i].set_occupancy(-1);
        // }

        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        std::vector<double> Iq(debye_axis.bins, 0);

        container::Container2D<double> distances(waters.size(), waters.size());
        for (unsigned int i = 0; i < waters.size(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                distances.index(i, j) = waters[i].distance(waters[j]);
                distances.index(j, i) = distances.index(i, j);
            }
        }

        auto ff = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::O);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            std::cout << std::round(double(q)/(q0+debye_axis.bins)*100) << "%\r" << std::flush;
            for (unsigned int i = 0; i < waters.size(); ++i) {
                for (unsigned int j = 0; j < waters.size(); ++j) {
                    double qd = constants::axes::q_vals[q]*distances.index(i, j);
                    if (qd < 1e-3) {
                        double qd2 = qd*qd;
                        Iq[q-q0] += waters[i].get_occupancy()*waters[j].get_occupancy()*(1 - qd2/6 + qd2*qd2/120);
                    } else {
                        Iq[q-q0] += waters[i].get_occupancy()*waters[j].get_occupancy()*std::sin(qd)/qd;
                    }
                }
            }
            Iq[q] *= std::pow(ff.evaluate(constants::axes::q_vals[q]), 2);
        }

        SimpleDataset dataset(debye_axis.as_vector(), Iq);
        dataset.save(settings::general::output + "waxsis_ww.dat");   
    }

    if (calc_density) {
        data::Molecule protein(env_ordered);
        Axis axis(0, 10, 50);

        auto calc_density = [&] (const data::Molecule& protein) {
            std::vector<double> min_dists;
            auto atoms = protein.get_atoms();
            for (auto& water : protein.get_waters()) {
                double shortest = std::numeric_limits<double>::max();

                for (unsigned int i = 0; i < protein.size_atom(); ++i) {
                    const auto& atom = atoms[i];
                    shortest = std::min(shortest, water.get_coordinates().distance2(atom.get_coordinates()));
                }
                min_dists.push_back((shortest < 0 ? -1 : 1) * std::sqrt(std::abs(shortest)));
            }

            hist::Histogram hist(axis);
            hist.bin(min_dists);
            hist.normalize_max();
            return hist;
        };

        hist::Histogram hist_radial, hist_axes, hist_jan;
        {
            settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
            data::Molecule tmp(env_ordered);
            tmp.generate_new_hydration();
            hist_radial = calc_density(tmp);
        }
        // {
        //     settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::AxesStrategy;
        //     data::Molecule tmp(env_ordered);
        //     tmp.generate_new_hydration();
        //     hist_axes = calc_density(tmp);
        // }
        // {
        //     settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::JanStrategy;
        //     data::Molecule tmp(env_ordered); 
        //     tmp.generate_new_hydration();
        //     hist_jan = calc_density(tmp);
        // }

        auto hist_ordered = calc_density(ordered);
        // auto hist_disordered = calc_density(disordered);

        plots::PlotDataset()
            .plot(hist_ordered.as_dataset(), plots::PlotOptions(style::draw::points, {{"xlabel", "d $[\\AA]$"}, {"ylabel", "Density [arb]"}, {"title", name}, {"legend", "WAXSiS"}, {"color", "k"}, {"lines", true}, {"lw", 0.5}}))
            // .plot(hist_disordered.as_dataset(), plots::PlotOptions(style::draw::points, {{"legend", "disordered"}, {"color", "r"}, {"lines", true}, {"lw", 0.5}}))
            .plot(hist_radial.as_dataset(), plots::PlotOptions(style::draw::points, {{"legend", "AUSAXS"}, {"color", "b"}, {"lines", true}, {"lw", 0.5}}))
            // .plot(hist_axes.as_dataset(), plots::PlotOptions(style::draw::points, {{"legend", "axes"}, {"color", "g"}, {"lines", true}, {"lw", 0.5}}))
            // .plot(hist_jan.as_dataset(), plots::PlotOptions(style::draw::points, {{"legend", "jan"}, {"color", "y"}, {"lines", true}, {"lw", 0.5}}))
        .save(settings::general::output + name + ".png");
    }

    return 0;
}