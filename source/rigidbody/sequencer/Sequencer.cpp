// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <io/ExistingFile.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hydrate/ExplicitHydration.h>
#include <data/atoms/AtomMetadata.h>
#include <utility/Console.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

Sequencer::Sequencer() : LoopElement(nullptr, 1), setup_loop(this), rigidbody(nullptr) {
    data::AtomMetadata::store_backbone = true;
    data::AtomMetadata::store_residue_seq = true;
}

Sequencer::Sequencer(const io::ExistingFile& saxs) : LoopElement(nullptr, 1), setup_loop(this, saxs), rigidbody(nullptr) {
    data::AtomMetadata::store_backbone = true;
    data::AtomMetadata::store_residue_seq = true;
}

Sequencer::~Sequencer() = default;

LoopElement& Sequencer::end() {
    throw ausaxs::except::runtime_error("Sequencer::end: Too many end() calls detected.");
}

void Sequencer::run() {
    throw ausaxs::except::logic_error("Sequencer::run: Use execute() to run the sequencer. Calling run() directly skips rigidbody and controller initialization.");
}

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

observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> Sequencer::_get_best_conf() const {
    assert(rigidbody != nullptr && "Sequencer::_get_best_conf: Rigidbody not set.");
    assert(_get_controller() != nullptr && "Sequencer::_get_best_conf: Controller not set.");
    return _get_controller()->get_current_best_config();
}

observer_ptr<controller::IController> Sequencer::_get_controller() const {
    assert(rigidbody != nullptr && "Sequencer::_get_controller: Rigidbody not set.");
    assert(rigidbody->controller != nullptr && "Sequencer::_get_controller: Controller not set.");
    return rigidbody->controller.get();
}

SetupElement& Sequencer::setup() {return setup_loop;}

std::shared_ptr<fitter::FitResult> Sequencer::execute() {
    _clear_stop_request(); // a stop requested while nothing was running must not immediately kill this run
    _reset_counters();     // a previous run must not leak into this one
    _recount_total_iterations(this);
    auto saxs_path = setup()._get_saxs_path();
    if (!saxs_path.exists()) {throw ausaxs::except::runtime_error("Sequencer::execute: SAXS file \"" + saxs_path.str() + "\" does not exist.");}
    rigidbody->molecule.generate_new_hydration(); // some setup elements requires access to the hydration generators

    // run the setup elements, defining all of the necessary parameters
    for (auto& e : setup()._get_elements()) {
        e->run();
    }

    // prepare the fitter for the actual optimization
    _get_controller()->setup(saxs_path);
    for (auto& e : LoopElement::elements) {
        if (_stop_requested()) {break;}
        e->run();
    }

    // a stopped run is still finished off normally below, so the caller gets the best conformation found so far
    if (_stop_requested()) {
        console::print_text_critical("Refinement stopped. Returning the best conformation found so far.");
    }

    // restore the best hydration shell before the final fit
    auto best_conf = _get_best_conf();
    if (!best_conf->waters.empty()) {
        rigidbody->molecule.clear_hydration();
        auto hydration = std::make_unique<hydrate::ExplicitHydration>(std::move(best_conf->waters));
        hydration->expanded_across_symmetry = true; // snapshotted from a molecule already hydrated by generate_new_hydration()
        rigidbody->molecule.get_body(0).set_hydration(std::move(hydration));
    }

    // update the fitter with the restored hydration shell and symmetry state. This is unconditional because the last
    // step may have been rejected, in which case finish_step() left the fitter holding the rejected candidate's model.
    _get_controller()->get_fitter()->set_model(rigidbody->molecule.get_histogram());

    return _get_controller()->get_fitter()->unconstrained_fit();
}