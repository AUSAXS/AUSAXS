#include <data/DataFwd.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace ausaxs::symmetry::detail {
    // Helper struct to store the compact coordinates for a body
    // Indexing is as follows: [symmetry #][replication #]
    // The first index, [0][0], is reserved for the original coordinates
    struct CompactCoordinateSymmetries {
        std::vector<std::vector<hist::detail::CompactCoordinates>> atomic;
        std::vector<std::vector<hist::detail::CompactCoordinates>> waters;
    };

    std::vector<CompactCoordinateSymmetries> generate_transformed_data(const data::Molecule& protein);
    CompactCoordinateSymmetries generate_transformed_data(const data::Body& body);
    std::vector<std::vector<hist::detail::CompactCoordinates>> generate_transformed_waters(const data::Body& body);
}