// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <io/ExistingFile.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

Sequencer::Sequencer() : LoopElement(nullptr, 1), setup_loop(this), rigidbody(nullptr) {}

Sequencer::Sequencer(const io::ExistingFile& saxs) : LoopElement(nullptr, 1), setup_loop(this, saxs), rigidbody(nullptr) {}

Sequencer::~Sequencer() = default;

observer_ptr<rigidbody::Rigidbody> Sequencer::_get_rigidbody() const {
    assert(rigidbody != nullptr && "Sequencer::_get_rigidbody: Rigidbody not set.");
    return rigidbody;
}

void Sequencer::_set_rigidbody(observer_ptr<Rigidbody> rigidbody) {
    assert(rigidbody != nullptr && "Sequencer::_set_rigidbody: Rigidbody must not be null.");
    setup_loop._set_active_body(rigidbody);
}

observer_ptr<data::Molecule> Sequencer::_get_molecule() const {
    assert(rigidbody != nullptr && "Sequencer::_get_molecule: Rigidbody not set.");
    return &rigidbody->molecule;
}

observer_ptr<const Sequencer> Sequencer::_get_sequencer() const {
    return this;
}

observer_ptr<Sequencer> Sequencer::_get_sequencer() {
    return this;
}

observer_ptr<rigidbody::detail::BestConf> Sequencer::_get_best_conf() const {
    assert(rigidbody != nullptr && "Sequencer::_get_best_conf: Rigidbody not set.");
    assert(_get_controller() != nullptr && "Sequencer::_get_best_conf: Controller not set.");
    return _get_controller()->current_best();
}

observer_ptr<controller::IController> Sequencer::_get_controller() const {
    assert(rigidbody != nullptr && "Sequencer::_get_controller: Rigidbody not set.");
    assert(rigidbody->controller != nullptr && "Sequencer::_get_controller: Controller not set.");
    return rigidbody->controller.get();
}

observer_ptr<SetupElement> Sequencer::setup() {return &setup_loop;}

std::shared_ptr<fitter::FitResult> Sequencer::execute() {
    auto saxs_path = setup()->_get_saxs_path();
    if (!saxs_path.exists()) {throw std::runtime_error("Sequencer::execute: SAXS file \"" + saxs_path.str() + "\" does not exist.");}
    rigidbody->molecule.generate_new_hydration(); // some setup elements requires access to the hydration generators

    // run the setup elements, defining all of the necessary parameters
    for (auto& e : setup()->_get_elements()) {
        e->run();
    }

    // prepare the fitter for the actual optimization
    _get_controller()->setup(saxs_path);
    for (auto& e : LoopElement::elements) {
        e->run();
    }

    return _get_controller()->get_fitter()->unconstrained_fit();
}