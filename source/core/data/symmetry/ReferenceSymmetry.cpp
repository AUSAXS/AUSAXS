// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/ReferenceSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <constants/Constants.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // The orientation every participant of a shared symmetry generates its copies in. It must be one common rotation rather than each body's own: the copies
    // are placed by a single operator acting on the whole asymmetric unit, so per-body operators would scatter the participants of a copy relative to each
    // other and the assembly would no longer be symmetric. The primary body owns the symmetry, so its orientation is the natural choice.
    std::optional<Matrix<double>> shared_orientation(observer_ptr<const data::Molecule> molecule, int primary_body) {
        return molecule->get_body(primary_body).symmetry().get_obj()->orientation;
    }
}

ReferenceSymmetry::ReferenceSymmetry(std::unique_ptr<ISymmetry> base, std::vector<int> bodies, std::vector<int> slots, observer_ptr<const data::Molecule> molecule)
    : base(std::move(base)), bodies(std::move(bodies)), slots(std::move(slots)), molecule(molecule)
{
    assert(this->base != nullptr && "ReferenceSymmetry: base symmetry cannot be null.");
    assert(!this->bodies.empty() && "ReferenceSymmetry: at least one participating body is required.");
    assert(this->bodies.size() == this->slots.size() && "ReferenceSymmetry: bodies and slots must be parallel.");
    assert(this->molecule != nullptr && "ReferenceSymmetry: a molecule is required to determine the combined centre of mass.");
}

Vector3<double> ReferenceSymmetry::combined_cm() const {
    Vector3<double> sum{0, 0, 0};
    double total = 0;
    for (int idx : bodies) {
        const auto& body = molecule->get_body(idx);
        double mass = 0;
        for (const auto& a : body.get_atoms()) {mass += constants::mass::get_mass(a.form_factor_type());}
        sum += body.get_cm()*mass;
        total += mass;
    }
    assert(0 < total && "ReferenceSymmetry::combined_cm: participating bodies contain no mass.");
    return sum/total;
}

AffineTransform ReferenceSymmetry::_make_transform(const Vector3<double>& anchor, int rep) const {
    return base->_get_transform(anchor, rep);
}

// The per-body cm is ignored: every participating body must rotate about the shared combined centre so the whole assembly is replicated as one rigid unit.
// Note that the centre is read off the participants' current positions, so it drifts whenever they move relative to one another rather than as a rigid group.
// That is deliberate: pinning it to the conformation the symmetry was defined in would couple the shared parameters to every real-space move of a participant,
// which is exactly what stops the two sets of parameters from being optimised at the same time. The copies remain mutually congruent through the drift.
Vector3<double> ReferenceSymmetry::_transform_anchor(const Vector3<double>&) const {return combined_cm();}

std::optional<Matrix<double>> ReferenceSymmetry::_transform_orientation(const std::optional<Matrix<double>>&) const {
    return shared_orientation(molecule, bodies.front());
}

unsigned int ReferenceSymmetry::repetitions() const {return base->repetitions();}
bool ReferenceSymmetry::is_closed() const {return base->is_closed();}
std::string ReferenceSymmetry::type_name() const {return base->type_name();}
std::span<double> ReferenceSymmetry::span_translation() {return base->span_translation();}
std::span<double> ReferenceSymmetry::span_rotation() {return base->span_rotation();}
std::vector<SymmetricDuplicatePair> ReferenceSymmetry::internal_pair_schedule() const {return base->internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetry::clone() const {
    return std::make_unique<ReferenceSymmetry>(base->clone(), bodies, slots, molecule);
}

ISymmetry& ReferenceSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const ReferenceSymmetry*>(other);
    assert(cast != nullptr && "Can only add ReferenceSymmetry with another ReferenceSymmetry.");
    base->add(cast->base.get());
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

// the view must place its copies with exactly the operator the target does, so it borrows the target's anchor rather than deriving one of its own
AffineTransform ReferenceSymmetryView::_make_transform(const Vector3<double>& anchor, int rep) const {
    return target()->_make_transform(anchor, rep);
}

Vector3<double> ReferenceSymmetryView::_transform_anchor(const Vector3<double>&) const {return target()->combined_cm();}

std::optional<Matrix<double>> ReferenceSymmetryView::_transform_orientation(const std::optional<Matrix<double>>&) const {
    return shared_orientation(molecule, primary_body);
}

unsigned int ReferenceSymmetryView::repetitions() const {return target()->repetitions();}
bool ReferenceSymmetryView::is_closed() const {return target()->is_closed();}
std::string ReferenceSymmetryView::type_name() const {return target()->type_name();}
std::vector<SymmetricDuplicatePair> ReferenceSymmetryView::internal_pair_schedule() const {return target()->internal_pair_schedule();}

std::unique_ptr<ISymmetry> ReferenceSymmetryView::clone() const {
    return std::make_unique<ReferenceSymmetryView>(molecule, primary_body, symmetry_index);
}

// the shared symmetry is perturbed once, through the primary body that owns it; the view is inert to the optimiser so the shared parameters are not updated twice
std::span<double> ReferenceSymmetryView::span_translation() {return {};}
std::span<double> ReferenceSymmetryView::span_rotation() {return {};}

ISymmetry& ReferenceSymmetryView::add(observer_ptr<const ISymmetry>) {return *this;}
