// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PairSchedule.h>

#include <cassert>
#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::symmetry;

CompositeSymmetry::CompositeSymmetry(std::unique_ptr<ISymmetry> inner, std::unique_ptr<ISymmetry> outer)
    : inner(std::move(inner)), outer(std::move(outer))
{
    assert(this->inner != nullptr && this->outer != nullptr && "CompositeSymmetry: sub-symmetries cannot be null.");
}

unsigned int CompositeSymmetry::repetitions() const {
    return (1+inner->repetitions())*(1+outer->repetitions()) - 1;
}

bool CompositeSymmetry::is_closed() const {
    return inner->is_closed() && outer->is_closed();
}

std::unique_ptr<ISymmetry> CompositeSymmetry::clone() const {
    return std::make_unique<CompositeSymmetry>(inner->clone(), outer->clone());
}

std::function<Vector3<double>(Vector3<double>)> CompositeSymmetry::get_transform(const Vector3<double>& cm, int rep) const {
    assert(0 < rep && rep <= static_cast<int>(repetitions()) && "CompositeSymmetry::get_transform: repetition index out of range.");

    // copy `rep` decodes to (outer copy k, inner copy j); the inner unit is replicated by the outer
    int stride = 1 + static_cast<int>(inner->repetitions());
    int k = rep / stride;
    int j = rep % stride;

    auto identity = [](Vector3<double> v) {return v;};
    std::function<Vector3<double>(Vector3<double>)> inner_t = identity;
    std::function<Vector3<double>(Vector3<double>)> outer_t = identity;
    if (j != 0) {inner_t = inner->get_transform(cm, j);}
    if (k != 0) {outer_t = outer->get_transform(cm, k);}

    return [inner_t = std::move(inner_t), outer_t = std::move(outer_t)](Vector3<double> v) {
        return outer_t(inner_t(v));
    };
}

std::vector<SymmetricDuplicatePair> CompositeSymmetry::internal_pair_schedule() const {
    // The partition of copy-pairs into equal-distance classes is invariant to the body cm:
    // both sub-symmetries share it, so the whole composite is covariant under conjugation by
    // Trans(cm), which preserves the relative transforms the bucketer keys on. Any fixed cm
    // therefore yields the correct partition; the offsets/angles of the sub-symmetries do
    // matter, so this is recomputed on each call rather than cached.
    Vector3<double> cm{0, 0, 0};
    int n = static_cast<int>(repetitions()) + 1;

    std::vector<AffineTransform> placements;
    placements.reserve(n);
    placements.push_back({Matrix<double>::identity(3), {0, 0, 0}}); // placement 0 = original body
    for (int p = 1; p < n; ++p) {
        // recover the affine map (R, T) by probing the transform: T = f(0), R*e_c = f(e_c) - T
        auto f = get_transform(cm, p);
        Vector3<double> T = f({0, 0, 0});
        Matrix<double> R(3, 3);
        for (int c = 0; c < 3; ++c) {
            Vector3<double> e{c == 0 ? 1.0 : 0.0, c == 1 ? 1.0 : 0.0, c == 2 ? 1.0 : 0.0};
            Vector3<double> col = f(e) - T;
            R(0, c) = col.x();
            R(1, c) = col.y();
            R(2, c) = col.z();
        }
        placements.push_back({std::move(R), std::move(T)});
    }
    return compute_pair_schedule(placements);
}

// a composite has two parameter sets; std::span cannot describe both at once. Callers must reach
// the sub-symmetries via for_each_leaf instead, so calling these directly is a programming error.
std::span<double> CompositeSymmetry::span_translation() {
    throw std::runtime_error("CompositeSymmetry::span_translation: a composite has no single contiguous parameter span; use symmetry::for_each_leaf to reach its sub-symmetries.");
}
std::span<double> CompositeSymmetry::span_rotation() {
    throw std::runtime_error("CompositeSymmetry::span_rotation: a composite has no single contiguous parameter span; use symmetry::for_each_leaf to reach its sub-symmetries.");
}

void ausaxs::symmetry::for_each_leaf(ISymmetry& sym, const std::function<void(ISymmetry&)>& fn) {
    if (auto* composite = dynamic_cast<CompositeSymmetry*>(&sym)) {
        for_each_leaf(*composite->inner, fn);
        for_each_leaf(*composite->outer, fn);
        return;
    }
    fn(sym);
}

ISymmetry& CompositeSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const CompositeSymmetry*>(other);
    assert(cast != nullptr && "Can only add CompositeSymmetry with another CompositeSymmetry.");
    inner->add(cast->inner.get());
    outer->add(cast->outer.get());
    return *this;
}
