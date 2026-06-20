// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/ReferenceSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

ReferenceSymmetry::ReferenceSymmetry(CyclicSymmetry base, std::vector<int> bodies, std::vector<int> slots, observer_ptr<const data::Molecule> molecule)
    : base(std::move(base)), bodies(std::move(bodies)), slots(std::move(slots)), molecule(molecule)
{
    assert(!this->bodies.empty() && "ReferenceSymmetry: at least one participating body is required.");
    assert(this->bodies.size() == this->slots.size() && "ReferenceSymmetry: bodies and slots must be parallel.");
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
    // the per-body cm is ignored: every participating body must rotate about the shared combined centre so the whole assembly is replicated as one rigid unit
    return base.get_transform(combined_cm(), rep);
}

unsigned int ReferenceSymmetry::repetitions() const {return base.repetitions();}
bool ReferenceSymmetry::is_closed() const {return base.is_closed();}
std::span<double> ReferenceSymmetry::span_translation() {return base.span_translation();}
std::span<double> ReferenceSymmetry::span_rotation() {return base.span_rotation();}
std::vector<SymmetricDuplicatePair> ReferenceSymmetry::internal_pair_schedule() const {return base.internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetry::clone() const {
    return std::make_unique<ReferenceSymmetry>(base, bodies, slots, molecule);
}

ISymmetry& ReferenceSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const ReferenceSymmetry*>(other);
    assert(cast != nullptr && "Can only add ReferenceSymmetry with another ReferenceSymmetry.");
    base.add(&cast->base);
    return *this;
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
std::vector<SymmetricDuplicatePair> ReferenceSymmetryView::internal_pair_schedule() const {return target()->internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetryView::clone() const {
    return std::make_unique<ReferenceSymmetryView>(molecule, primary_body, symmetry_index);
}

// the shared symmetry is perturbed once, through the primary body that owns it; the view is inert to the optimiser so the shared parameters are not updated twice
std::span<double> ReferenceSymmetryView::span_translation() {return {};}
std::span<double> ReferenceSymmetryView::span_rotation() {return {};}

ISymmetry& ReferenceSymmetryView::add(observer_ptr<const ISymmetry>) {return *this;}
