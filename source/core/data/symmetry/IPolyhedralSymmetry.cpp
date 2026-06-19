// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/IPolyhedralSymmetry.h>
#include <data/symmetry/PairSchedule.h>
#include <math/MatrixUtils.h>

#include <array>
#include <cassert>
#include <cmath>
#include <set>
#include <typeinfo>

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
}

IPolyhedralSymmetry::GroupData IPolyhedralSymmetry::build(const std::vector<Matrix<double>>& generators, int order) {
    auto elements = close_group(generators, order);

    // distance-reuse schedule: the equivalence classes depend only on the fixed group structure,
    // not on the optimisable offset/frame, so the placements are the bare group rotations about the origin.
    std::vector<AffineTransform> placements;
    placements.reserve(elements.size());
    for (const auto& M : elements) {placements.push_back({M, {0, 0, 0}});}
    return GroupData{std::move(elements), compute_pair_schedule(placements)};
}

unsigned int IPolyhedralSymmetry::repetitions() const {
    return static_cast<unsigned int>(group().elements.size()) - 1;
}

bool IPolyhedralSymmetry::is_closed() const {return false;}

std::function<Vector3<double>(Vector3<double>)> IPolyhedralSymmetry::get_transform(const Vector3<double>& cm, int rep) const {
    const auto& G = group().elements;
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

std::span<double> IPolyhedralSymmetry::span_translation() {return std::span<double>(translation.begin(), translation.end());}
std::span<double> IPolyhedralSymmetry::span_rotation() {return std::span<double>(rotation.begin(), rotation.end());}

std::vector<SymmetricDuplicatePair> IPolyhedralSymmetry::internal_pair_schedule() const {
    return group().schedule;
}

ISymmetry& IPolyhedralSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const IPolyhedralSymmetry*>(other);
    assert(cast != nullptr && "Can only add PolyhedralSymmetry with another PolyhedralSymmetry.");
    assert(typeid(*this) == typeid(*cast) && "Cannot add PolyhedralSymmetry objects of different groups.");
    this->translation += cast->translation;
    this->rotation += cast->rotation;
    return *this;
}