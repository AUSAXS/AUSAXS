// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/RigidBody.h>
#include <fitter/SmartFitter.h>
#include <grid/Grid.h>
#include <io/ExistingFile.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

Sequencer::Sequencer() : LoopElement(nullptr, 1), setup_loop(this), rigidbody(nullptr), best(nullptr) {}

Sequencer::Sequencer(const io::ExistingFile& saxs) : LoopElement(nullptr, 1), setup_loop(this, saxs), rigidbody(nullptr), best(nullptr) {}

Sequencer::~Sequencer() = default;

observer_ptr<rigidbody::RigidBody>& Sequencer::_get_rigidbody() {
    return rigidbody;
}

observer_ptr<rigidbody::RigidBody> Sequencer::_get_rigidbody() const {
    return rigidbody;
}

observer_ptr<rigidbody::detail::BestConf> Sequencer::_get_best_conf() const {
    return best.get();
}

observer_ptr<const Sequencer> Sequencer::_get_sequencer() const {
    return this;
}

bool Sequencer::_optimize_step() const {
    return rigidbody->optimize_step(*best);
}

observer_ptr<SetupElement> Sequencer::setup() {return &setup_loop;}

std::shared_ptr<fitter::FitResult> Sequencer::execute() {
    auto saxs_path = setup()->_get_saxs_path();
    if (!saxs_path.exists()) {throw std::runtime_error("Sequencer::execute: SAXS file \"" + saxs_path.str() + "\"does not exist.");}
    rigidbody->generate_new_hydration(); // some setup elements requires access to the hydration generators

    // run the setup elements, defining all of the necessary parameters
    for (auto& e : setup()->_get_elements()) {
        e->run();
    }

    // prepare the fitter for the actual optimization
    rigidbody->prepare_fitter(saxs_path);
    best = std::make_unique<detail::BestConf>(std::make_shared<grid::Grid>(*rigidbody->get_grid()), rigidbody->get_waters(), rigidbody->fitter->fit_chi2_only());

    for (auto& e : LoopElement::elements) {
        e->run();
    }

    return rigidbody->get_unconstrained_fitter(saxs_path)->fit();
}