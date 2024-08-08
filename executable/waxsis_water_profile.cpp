#include <CLI/CLI.hpp>

#include <grid/Grid.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/All.h>
#include <constants/Constants.h>
#include <table/ArrayDebyeTable.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <dataset/SimpleDataset.h>
#include <form_factor/FormFactor.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);

    io::Folder input;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("folder", input, "Path to the MD SAXS folder.")->required()->check(CLI::ExistingDirectory);
    CLI11_PARSE(app, argc, argv);

    settings::grid::scaling = 1;
    settings::grid::width = 0.5;
    settings::grid::rvol = 2.5; // water + atom
    settings::molecule::use_effective_charge = false;
    settings::general::output = "output/";

    io::ExistingFile env_unordered(input + "excludedvolume_0.pdb");
    io::ExistingFile env_ordered(input + "prot+solvlayer_0.pdb");

    data::Molecule disordered(env_unordered);
    data::Molecule ordered(env_ordered);
    auto ordered_layer = ordered.get_waters();

    ordered.clear_hydration();
    auto grid = ordered.get_grid();
    std::vector<data::record::Water> disordered_waters;
    for (auto& water : disordered.get_waters()) {
        auto bin = grid->to_bins(water.get_coordinates());
        if (grid->grid.is_empty(bin.x(), bin.y(), bin.z())) {
            disordered_waters.emplace_back(water);
        }
    }
    std::vector<data::record::Water> ordered_waters;
    for (auto& water : ordered_layer) {
        auto bin = grid->to_bins(water.get_coordinates());
        if (grid->grid.is_empty(bin.x(), bin.y(), bin.z())) {
            ordered_waters.emplace_back(water);
        }
    }

    ordered = data::Molecule(std::vector<data::record::Atom>{}, ordered_waters);
    disordered = data::Molecule(std::vector<data::record::Atom>{}, disordered_waters);

    ordered.save(settings::general::output + "ordered.pdb");
    disordered.save(settings::general::output + "disordered.pdb");

    //#######################//
    //### CALC SCATTERING ###//
    //#######################//
    auto sinqd_table = table::ArrayDebyeTable::get_default_table();
    std::vector<data::record::Water> waters = ordered.get_waters();
    waters.insert(waters.end(), disordered.get_waters().begin(), disordered.get_waters().end());
    for (unsigned int i = disordered.size_water(); i < waters.size(); ++i) {
        waters[i].set_occupancy(-1);
    }

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
    return 0;
}