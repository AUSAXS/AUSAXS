// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PolyhedralSymmetry.h>
#include <data/symmetry/PairSchedule.h>
#include <math/MatrixUtils.h>

#include <array>
#include <cassert>
#include <cmath>
#include <numbers>
#include <set>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // hashable key for a rotation matrix, used to deduplicate group elements
    std::array<long, 9> matrix_key(const Matrix<double>& M) {
        std::array<long, 9> k;
        int idx = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                k[idx++] = std::llround(M(i, j)*1e6);
            }
        }
        return k;
    }

    // generate a finite rotation group as the closure of a set of generators (BFS).
    // element 0 is always the identity; the closure is hard-capped so inconsistent
    // generators fail the assert below rather than looping forever.
    std::vector<Matrix<double>> close_group(const std::vector<Matrix<double>>& generators, int expected_order) {
        std::vector<Matrix<double>> elements;
        std::set<std::array<long, 9>> seen;
        auto add = [&](const Matrix<double>& M) {
            if (seen.insert(matrix_key(M)).second) {elements.push_back(M);}
        };
        add(Matrix<double>::identity(3));
        for (const auto& g : generators) {add(g);}
        for (std::size_t i = 0; i < elements.size() && static_cast<int>(elements.size()) <= expected_order; ++i) {
            for (const auto& g : generators) {add(elements[i]*g);}
        }
        assert(static_cast<int>(elements.size()) == expected_order && "PolyhedralSymmetry: group closure produced an unexpected order");
        return elements;
    }

    std::vector<Matrix<double>> build_group(PolyhedralGroup group) {
        constexpr double pi = std::numbers::pi;
        auto rot = [](Vector3<double> axis, double angle) {return matrix::rotation_matrix<double>(axis, angle);};
        switch (group) {
            // 3-fold body-diagonal rotation + 2-fold face rotation generate the rotation group A4
            case PolyhedralGroup::tetrahedral:
                return close_group({rot({1, 1, 1}, 2*pi/3), rot({0, 0, 1}, pi)}, 12);
            // 4-fold face rotation + 3-fold body-diagonal rotation generate the rotation group S4
            case PolyhedralGroup::octahedral:
                return close_group({rot({0, 0, 1}, pi/2), rot({1, 1, 1}, 2*pi/3)}, 24);
            // 5-fold vertex + 2-fold edge + 3-fold face rotations of an icosahedron with
            // vertices at the cyclic permutations of (0, +-1, +-phi) generate the group A5
            case PolyhedralGroup::icosahedral: {
                double phi = (1 + std::sqrt(5.0))/2;
                return close_group({rot({0, 1, phi}, 2*pi/5), rot({1, 0, 0}, pi), rot({1, 1, 1}, 2*pi/3)}, 60);
            }
        }
        assert(false && "PolyhedralSymmetry: unknown group");
        return {};
    }

    // group rotation matrices, built once per process; element 0 is the identity
    const std::vector<Matrix<double>>& group_elements(PolyhedralGroup group) {
        static const std::vector<Matrix<double>> t = build_group(PolyhedralGroup::tetrahedral);
        static const std::vector<Matrix<double>> o = build_group(PolyhedralGroup::octahedral);
        static const std::vector<Matrix<double>> i = build_group(PolyhedralGroup::icosahedral);
        switch (group) {
            case PolyhedralGroup::tetrahedral: return t;
            case PolyhedralGroup::octahedral:  return o;
            case PolyhedralGroup::icosahedral: return i;
        }
        assert(false && "PolyhedralSymmetry: unknown group");
        return t;
    }

    // distance-reuse schedule, built once per process. The equivalence classes depend only
    // on the fixed group structure, not on the optimisable offset/frame, so the placements
    // are taken as the bare group rotations about the origin.
    const std::vector<CopyPair>& group_schedule(PolyhedralGroup group) {
        auto build = [](PolyhedralGroup g) {
            const auto& G = group_elements(g);
            std::vector<AffineTransform> placements;
            placements.reserve(G.size());
            for (const auto& M : G) {placements.push_back({M, {0, 0, 0}});}
            return compute_pair_schedule(placements);
        };
        static const std::vector<CopyPair> t = build(PolyhedralGroup::tetrahedral);
        static const std::vector<CopyPair> o = build(PolyhedralGroup::octahedral);
        static const std::vector<CopyPair> i = build(PolyhedralGroup::icosahedral);
        switch (group) {
            case PolyhedralGroup::tetrahedral: return t;
            case PolyhedralGroup::octahedral:  return o;
            case PolyhedralGroup::icosahedral: return i;
        }
        assert(false && "PolyhedralSymmetry: unknown group");
        return t;
    }
}

PolyhedralSymmetry::PolyhedralSymmetry(PolyhedralGroup group) : group(group) {}

unsigned int PolyhedralSymmetry::repetitions() const {
    return static_cast<unsigned int>(group_elements(group).size()) - 1;
}

bool PolyhedralSymmetry::is_closed() const {return false;}

std::unique_ptr<ISymmetry> PolyhedralSymmetry::clone() const {
    return std::make_unique<PolyhedralSymmetry>(*this);
}

std::function<Vector3<double>(Vector3<double>)> PolyhedralSymmetry::get_transform(const Vector3<double>& cm, int rep) const {
    const auto& G = group_elements(group);
    assert(0 < rep && rep < static_cast<int>(G.size()) && "PolyhedralSymmetry::get_transform: repetition index out of range.");

    // copy `rep` is  v -> c + F G_rep F^T (v - c),  with c = cm + offset and F the frame orientation
    Matrix<double> F = matrix::rotation_matrix<double>(rotation);
    Matrix<double> R = F*G[rep]*F.transpose();
    Vector3<double> c = cm + translation;
    Vector3<double> T = c - R*c;
    return [R = std::move(R), T = std::move(T)](Vector3<double> v) {
        return R*v + T;
    };
}

std::span<double> PolyhedralSymmetry::span_translation() {return std::span<double>(translation.begin(), translation.end());}
std::span<double> PolyhedralSymmetry::span_rotation() {return std::span<double>(rotation.begin(), rotation.end());}

std::vector<CopyPair> PolyhedralSymmetry::internal_pair_schedule() const {
    return group_schedule(group);
}

ISymmetry& PolyhedralSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const PolyhedralSymmetry*>(other);
    assert(cast != nullptr && "Can only add PolyhedralSymmetry with another PolyhedralSymmetry.");
    assert(cast->group == group && "Cannot add PolyhedralSymmetry objects of different groups.");
    this->translation += cast->translation;
    this->rotation += cast->rotation;
    return *this;
}
