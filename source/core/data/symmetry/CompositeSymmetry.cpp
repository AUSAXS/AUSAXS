// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/PairSchedule.h>

#include <cassert>
#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::symmetry;

CompositeSymmetry::CompositeSymmetry(std::unique_ptr<ISymmetry> inner, std::unique_ptr<ISymmetry> outer)
    : inner(std::move(inner)), outer(std::move(outer))
{
    assert(this->inner != nullptr && this->outer != nullptr && "CompositeSymmetry: sub-symmetries cannot be null.");

    // A composite anchors the whole nesting on the owning body's centre of mass, while a shared symmetry anchors on the combined centre of its participants.
    // Nesting one inside the other would compose two transforms built around different pivots, so the composite is restricted to per-body sub-symmetries.
    assert([&] {
        for (const auto* s : {this->inner.get(), this->outer.get()}) {
            if (dynamic_cast<const ReferenceSymmetry*>(s) || dynamic_cast<const ReferenceSymmetryView*>(s)) {return false;}
        }
        return true;
    }() && "CompositeSymmetry: a shared (reference) symmetry cannot be nested inside a composite.");
}

unsigned int CompositeSymmetry::repetitions() const {
    return (1+inner->repetitions())*(1+outer->repetitions()) - 1;
}

bool CompositeSymmetry::is_closed() const {
    return inner->is_closed() && outer->is_closed();
}

std::string CompositeSymmetry::type_name() const {
    return inner->type_name() + "-" + outer->type_name();
}

std::unique_ptr<ISymmetry> CompositeSymmetry::clone() const {
    return std::make_unique<CompositeSymmetry>(inner->clone(), outer->clone());
}

AffineTransform CompositeSymmetry::_make_transform(const Vector3<double>& cm, int rep) const {
    assert(0 < rep && rep <= static_cast<int>(repetitions()) && "CompositeSymmetry::_make_transform: repetition index out of range.");

    // copy `rep` decodes to (outer copy k, inner copy j); the inner unit is replicated by the outer
    int stride = 1 + static_cast<int>(inner->repetitions());
    int k = rep / stride;
    int j = rep % stride;

    AffineTransform inner_t, outer_t; // default-constructed to the identity
    if (j != 0) {inner_t = inner->_get_transform(cm, j);}
    if (k != 0) {outer_t = outer->_get_transform(cm, k);}

    // outer after inner: v -> R_o*(R_i*v + T_i) + T_o
    return {outer_t.rotation*inner_t.rotation, outer_t.rotation*inner_t.translation + outer_t.translation};
}

std::vector<SymmetricDuplicatePair> CompositeSymmetry::internal_pair_schedule() const {
    // The partition of copy-pairs into equal-distance classes is invariant to the body cm: both sub-symmetries share it, so the whole composite is covariant 
    // under conjugation by Trans(cm), which preserves the relative transforms the bucketer keys on. Any fixed cm therefore yields the correct partition; the 
    // offsets/angles of the sub-symmetries do matter, so this is recomputed on each call rather than cached.
    Vector3<double> cm{0, 0, 0};
    int n = static_cast<int>(repetitions()) + 1;

    std::vector<AffineTransform> placements;
    placements.reserve(n);
    placements.emplace_back(); // placement 0 = original body
    for (int p = 1; p < n; ++p) {placements.push_back(_get_transform(cm, p));}
    return compute_pair_schedule(placements);
}

// a composite has two parameter sets; std::span cannot describe both at once. Callers must reach the sub-symmetries via for_each_leaf instead, 
// so calling these directly is a programming error.
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
    // a ReferenceSymmetry is just a wrapper around its base (possibly itself composite); descend into it rather than treating the wrapper as an opaque leaf, 
    // so callers reach the real parameters regardless of what kind of symmetry is shared
    if (auto* ref = dynamic_cast<ReferenceSymmetry*>(&sym)) {
        for_each_leaf(*ref->base, fn);
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
