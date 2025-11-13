// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace ausaxs::symmetry::detail {
    // Helper struct to store the compact coordinates for a body
    // Indexing is as follows: [symmetry #][replication #]
    // The first index, [0][0], is reserved for the original coordinates
    template<bool variable_bin_width>
    struct BodySymmetryData {
        template<typename T>
        using symmetry_t = std::vector<T>;

        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the symmetries, the inner loop is over the repetitions
        // index [0][0] is the original data
        symmetry_t<repetition_t<hist::detail::CompactCoordinates<variable_bin_width>>> atomic;
    };

    template<bool variable_bin_width>
    struct SymmetryData {
        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the repetitions
        repetition_t<hist::detail::CompactCoordinates<variable_bin_width>> data;
    };

    template<bool variable_bin_width>
    std::pair<std::vector<BodySymmetryData<variable_bin_width>>, hist::detail::CompactCoordinates<variable_bin_width>> generate_transformed_data(const data::Molecule& protein);
    
    template<bool variable_bin_width>
    BodySymmetryData<variable_bin_width> generate_transformed_data(const data::Body& body);

    template<bool variable_bin_width>
    SymmetryData<variable_bin_width> generate_transformed_data(const data::Body& body, int isym);
}