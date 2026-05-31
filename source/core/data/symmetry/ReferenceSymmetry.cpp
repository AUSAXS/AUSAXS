// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/ReferenceSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/state/Signaller.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

ReferenceSymmetry::ReferenceSymmetry(CyclicSymmetry base, std::vector<int> bodies, observer_ptr<const data::Molecule> molecule)
    : base(std::move(base)), bodies(std::move(bodies)), molecule(molecule)
{
    assert(!this->bodies.empty() && "ReferenceSymmetry: at least one participating body is required.");
    assert(this->molecule != nullptr && "ReferenceSymmetry: a molecule is required to determine the combined centre of mass.");
}

Vector3<double> ReferenceSymmetry::combined_cm() const {
    // atom-count weighted average of the participating bodies' centres of mass
    Vector3<double> sum{0, 0, 0};
    std::size_t total = 0;
    for (int idx : bodies) {
        const auto& body = molecule->get_body(idx);
        std::size_t n = body.size_atom();
        sum += body.get_cm()*static_cast<double>(n);
        total += n;
    }
    assert(0 < total && "ReferenceSymmetry::combined_cm: participating bodies contain no atoms.");
    return sum/static_cast<double>(total);
}

std::function<Vector3<double>(Vector3<double>)> ReferenceSymmetry::get_transform(const Vector3<double>&, int rep) const {
    // the per-body cm is ignored: every participating body must rotate about the shared
    // combined centre so the whole assembly is replicated as one rigid unit
    return base.get_transform(combined_cm(), rep);
}

unsigned int ReferenceSymmetry::repetitions() const {return base.repetitions();}
bool ReferenceSymmetry::is_closed() const {return base.is_closed();}
std::span<double> ReferenceSymmetry::span_translation() {return base.span_translation();}
std::span<double> ReferenceSymmetry::span_rotation() {return base.span_rotation();}
std::vector<CopyPair> ReferenceSymmetry::internal_pair_schedule() const {return base.internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetry::clone() const {
    return std::make_unique<ReferenceSymmetry>(base, bodies, molecule);
}

ISymmetry& ReferenceSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const ReferenceSymmetry*>(other);
    assert(cast != nullptr && "Can only add ReferenceSymmetry with another ReferenceSymmetry.");
    base.add(&cast->base);
    return *this;
}

void ReferenceSymmetry::signal_modified(observer_ptr<const signaller::Signaller> host, int index) const {
    // flag the primary body that owns this symmetry...
    host->modified_symmetry(index);

    // ...and every other participating body, whose view delegates to this shared symmetry, so the
    // partial histogram manager recomputes their copies too. host is the primary's signaller and
    // index its slot, which is exactly what the views point back to.
    int primary = bodies.front();
    for (int b : bodies) {
        if (b == primary) {continue;}
        const auto& body = molecule->get_body(b);
        for (int j = 0; j < static_cast<int>(body.size_symmetry()); ++j) {
            auto view = dynamic_cast<const ReferenceSymmetryView*>(body.symmetry().get(j));
            if (view != nullptr && view->primary_body == primary && view->symmetry_index == index) {
                body.get_signaller()->modified_symmetry(j);
            }
        }
    }
}

ReferenceSymmetryView::ReferenceSymmetryView(observer_ptr<const data::Molecule> molecule, int primary_body, int symmetry_index)
    : molecule(molecule), primary_body(primary_body), symmetry_index(symmetry_index)
{
    assert(molecule != nullptr && "ReferenceSymmetryView: molecule cannot be null.");
    assert(0 <= primary_body && "ReferenceSymmetryView: primary body index must be non-negative.");
    assert(0 <= symmetry_index && "ReferenceSymmetryView: symmetry index must be non-negative.");
}

observer_ptr<const ReferenceSymmetry> ReferenceSymmetryView::target() const {
    auto sym = molecule->get_body(primary_body).symmetry().get(symmetry_index);
    auto ref = dynamic_cast<const ReferenceSymmetry*>(sym);
    assert(ref != nullptr && "ReferenceSymmetryView::target: the referenced symmetry is not a ReferenceSymmetry.");
    return ref;
}

std::function<Vector3<double>(Vector3<double>)> ReferenceSymmetryView::get_transform(const Vector3<double>& cm, int rep) const {
    return target()->get_transform(cm, rep);
}

unsigned int ReferenceSymmetryView::repetitions() const {return target()->repetitions();}
bool ReferenceSymmetryView::is_closed() const {return target()->is_closed();}
std::vector<CopyPair> ReferenceSymmetryView::internal_pair_schedule() const {return target()->internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetryView::clone() const {
    return std::make_unique<ReferenceSymmetryView>(molecule, primary_body, symmetry_index);
}

// the shared symmetry is perturbed once, through the primary body that owns it; the view is
// inert to the optimiser so the shared parameters are not updated twice
std::span<double> ReferenceSymmetryView::span_translation() {return {};}
std::span<double> ReferenceSymmetryView::span_rotation() {return {};}

ISymmetry& ReferenceSymmetryView::add(observer_ptr<const ISymmetry>) {return *this;}

void ReferenceSymmetryView::signal_modified(observer_ptr<const signaller::Signaller> host, int index) const {
    // flag this body, then cascade through the owning symmetry so the whole group is recomputed
    host->modified_symmetry(index);
    target()->signal_modified(molecule->get_body(primary_body).get_signaller().get(), symmetry_index);
}
