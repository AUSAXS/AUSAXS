// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>

#include <type_traits>
#include <vector>

namespace ausaxs::symmetry::detail {
    /**
     * @brief The atoms of one copy of a body: a single set, or with @a form_factors one set per active form factor type.
     */
    template<bool form_factors>
    using AtomicCoordinates = std::conditional_t<form_factors, std::vector<hist::detail::CompactCoordinates>, hist::detail::CompactCoordinates>;

    // Helper struct to store the compact coordinates for a body
    // Indexing is as follows: [symmetry #][replication #]
    // The first index, [0][0], is reserved for the original coordinates
    template<bool form_factors>
    struct BodySymmetryData {
        template<typename T>
        using symmetry_t = std::vector<T>;

        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the symmetries, the inner loop is over the repetitions
        // index [0][0] is the original data
        symmetry_t<repetition_t<AtomicCoordinates<form_factors>>> atomic;
    };

    template<bool form_factors>
    struct SymmetryData {
        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the repetitions
        repetition_t<AtomicCoordinates<form_factors>> data;
    };

    template<bool form_factors>
    std::pair<std::vector<BodySymmetryData<form_factors>>, hist::detail::CompactCoordinates> generate_transformed_data(const data::Molecule& protein);

    template<bool form_factors>
    BodySymmetryData<form_factors> generate_transformed_data(const data::Body& body);

    template<bool form_factors>
    SymmetryData<form_factors> generate_transformed_data(const data::Body& body, int isym);
}