// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <io/ExistingFile.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

Sequencer::Sequencer() : LoopElement(nullptr, 1), setup_loop(this), rigidbody(nullptr), best(nullptr) {}

Sequencer::Sequencer(const io::ExistingFile& saxs) : LoopElement(nullptr, 1), setup_loop(this, saxs), rigidbody(nullptr), best(nullptr) {}

Sequencer::~Sequencer() = default;

observer_ptr<rigidbody::Rigidbody>& Sequencer::_get_rigidbody() {
    return rigidbody;
}

observer_ptr<rigidbody::Rigidbody> Sequencer::_get_rigidbody() const {
    return rigidbody;
}

observer_ptr<data::Molecule> Sequencer::_get_molecule() const {
    return &rigidbody->molecule;
}

observer_ptr<rigidbody::detail::BestConf> Sequencer::_get_best_conf() const {
    return rigidbody->controller->current_best();
}

observer_ptr<const Sequencer> Sequencer::_get_sequencer() const {
    return this;
}

bool Sequencer::_optimize_step() const {
    return rigidbody->controller->run_step();
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
    rigidbody->controller->setup(saxs_path);
    best = std::make_unique<detail::BestConf>(
        std::make_shared<grid::Grid>(*rigidbody->molecule.get_grid()), 
        rigidbody->molecule.get_waters(), 
        rigidbody->controller->get_fitter()->fit_chi2_only()
    );

    for (auto& e : LoopElement::elements) {
        e->run();
    }

    return rigidbody->controller->get_fitter()->unconstrained_fit();
}