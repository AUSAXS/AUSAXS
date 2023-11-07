#include <hist/symmetry/PlaneSymmetry.h>
#include <hist/detail/CompactCoordinates.h>

using namespace hist::symmetry;

Vector3<double> PlaneSymmetry::mirror(const Vector3<double>& v) const {
    auto& n = plane.normal;
    return v - 2*(v.dot(n))*n;
}

std::vector<double> PlaneSymmetry::calculate_cross_terms(hist::detail::CompactCoordinates&, hist::detail::CompactCoordinates&) {
    // std::vector<double> cross_terms;
    // cross_terms.reserve(body1.get_size()*body2.get_size());

    // for (unsigned int i = 0; i < body1.get_size(); ++i) {
    //     for (unsigned int j = 0; j < body2.get_size(); ++j) {
    //         auto& v1 = body1[i];
    //         auto& v2 = body2[j];
    //         auto& n = plane.normal;

    //         // mirror v2 across the plane
    //         auto v2_mirrored = v2 - 2*(v2.dot(n))*n;

    //         // calculate the distance between v1 and v2_mirrored
    //         float dx = v1.x - v2_mirrored.x;
    //         float dy = v1.y - v2_mirrored.y;
    //         float dz = v1.z - v2_mirrored.z;
    //         float dist = std::sqrt(dx*dx + dy*dy + dz*dz);

    //         // calculate the weight
    //         float weight = v1.w*v2.w;

    //         // add the cross term
    //         cross_terms.push_back(weight*dist);
    //     }
    // }

    // return cross_terms;
    return {};
}