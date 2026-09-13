// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/PartialBinEstimate.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/detail/BinEstimate.h>

#include <array>
#include <cmath>

using namespace ausaxs;

template<bool variable_bin_width>
int hist::detail::required_partial_bin_count(const data::Molecule& protein) {
    std::vector<Vector3<double>> copies;
    for (const auto& body : protein.get_bodies()) {
        if (body.size_symmetry() == 0 || body.get_atoms().empty()) {continue;}

        Vector3<double> lo{
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        };
        Vector3<double> hi{
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (const auto& atom : body.get_atoms()) {
            const auto& coordinates = atom.coordinates();
            lo.x() = std::min(lo.x(), static_cast<double>(coordinates.x()));
            lo.y() = std::min(lo.y(), static_cast<double>(coordinates.y()));
            lo.z() = std::min(lo.z(), static_cast<double>(coordinates.z()));
            hi.x() = std::max(hi.x(), static_cast<double>(coordinates.x()));
            hi.y() = std::max(hi.y(), static_cast<double>(coordinates.y()));
            hi.z() = std::max(hi.z(), static_cast<double>(coordinates.z()));
        }

        std::array<Vector3<double>, 8> corners = {
            Vector3<double>{lo.x(), lo.y(), lo.z()}, Vector3<double>{lo.x(), lo.y(), hi.z()},
            Vector3<double>{lo.x(), hi.y(), lo.z()}, Vector3<double>{lo.x(), hi.y(), hi.z()},
            Vector3<double>{hi.x(), lo.y(), lo.z()}, Vector3<double>{hi.x(), lo.y(), hi.z()},
            Vector3<double>{hi.x(), hi.y(), lo.z()}, Vector3<double>{hi.x(), hi.y(), hi.z()}
        };

        const auto centre = body.get_cm();
        for (int isym = 0; isym < body.size_symmetry(); ++isym) {
            const auto* const symmetry = body.symmetry().get(isym);
            for (int repetition = 1; repetition <= symmetry->repetitions(); ++repetition) {
                const auto transform = body.symmetry().get_transform(isym, centre, repetition);
                for (const auto& corner : corners) {copies.push_back(transform(corner));}
            }
        }
    }
    return required_bin_count<variable_bin_width>(protein.iterate_atoms(), protein.iterate_waters(), copies);
}

int hist::detail::grown_partial_bin_count(int required) {
    constexpr double growth_margin = 0.1;
    constexpr int headroom = 2;
    return std::max(static_cast<int>(std::ceil(required*(1+growth_margin))), required+headroom);
}

template int hist::detail::required_partial_bin_count<true>(const data::Molecule&);
template int hist::detail::required_partial_bin_count<false>(const data::Molecule&);