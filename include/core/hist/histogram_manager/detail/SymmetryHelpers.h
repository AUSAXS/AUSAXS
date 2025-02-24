#pragma once

#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace ausaxs::symmetry::detail {
    // Helper struct to store the compact coordinates for a body
    // Indexing is as follows: [symmetry #][replication #]
    // The first index, [0][0], is reserved for the original coordinates
    struct BodySymmetryData {
        template<typename T>
        using symmetry_t = std::vector<T>;

        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the symmetries, the inner loop is over the repetitions
        // index [0][0] is the original data
        symmetry_t<repetition_t<hist::detail::CompactCoordinates>> atomic;
    };

    struct SymmetryData {
        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the repetitions
        repetition_t<hist::detail::CompactCoordinates> data;
    };

    std::pair<std::vector<BodySymmetryData>, hist::detail::CompactCoordinates> generate_transformed_data(const data::Molecule& protein);
    BodySymmetryData generate_transformed_data(const data::Body& body);
    SymmetryData generate_transformed_data(const data::Body& body, int isym);
}