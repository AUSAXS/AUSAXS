// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/culling/CounterCulling.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;

void hydrate::BodyCounterCulling::set_body_ratios(const std::vector<double>& body_ratios) {
    this->body_ratios = body_ratios;
}

void hydrate::BodyCounterCulling::cull(std::span<grid::GridMember<data::Water>>& placed_water) const {
    if (body_ratios.empty()) {return CounterCulling(molecule).cull(placed_water);}
    double total_reduction_factor = static_cast<double>(placed_water.size())/target_count;
    double total_weight = std::accumulate(body_ratios.begin(), body_ratios.end(), 0.0);
    assert(total_weight != 0 && "BodyCounterCulling::cull: Division by zero. The molecule has no atoms.");

    std::vector<Vector3<double>> cms(molecule->size_body());
    std::vector<double> reduction_factors(cms.size());
    for (unsigned int i = 0; i < molecule->size_body(); i++) {
        cms[i] = molecule->get_body(i).get_cm();
        reduction_factors[i] = total_reduction_factor*body_ratios[i]/total_weight;
    }

    std::vector<bool> removed_water(placed_water.size(), false);
    std::vector<unsigned int> body_counts(molecule->size_body(), 0);
    std::vector<double> body_counters(molecule->size_body(), 0);
    unsigned int rm_index = 0; // current index in removed_water
    unsigned int fw_index = 0; // current index in final_waters
    for (const auto& water : placed_water) {
        const Vector3<double>& pos = water.get_absolute_loc();
        double min_dist = std::numeric_limits<double>::max();
        unsigned int min_index = 0;
        for (unsigned int i = 0; i < cms.size(); ++i) {
            double dist = pos.distance2(cms[i]);
            if (dist < min_dist) {
                min_dist = dist;
                min_index = i;
            }
        }

        body_counters[min_index] += reduction_factors[min_index];
        if (body_counts[min_index] < body_counters[min_index]) {
            body_counts[min_index]++;
            removed_water[fw_index++] = true;
            continue;
        }
        removed_water[rm_index++] = true;
    }
    molecule->get_grid()->remove_waters(removed_water);
}