/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/RigidBody.h>
#include <fitter/LinearFitter.h>
#include <data/record/Water.h>
#include <hydrate/Grid.h>
#include <io/ExistingFile.h>

using namespace rigidbody::sequencer;

Sequencer::Sequencer(const io::ExistingFile& saxs, observer_ptr<RigidBody> rigidbody) : LoopElement(this, 1), rigidbody(rigidbody) {
    rigidbody->generate_new_hydration();
    rigidbody->prepare_fitter(saxs);
    best = std::make_unique<detail::BestConf>(std::make_shared<grid::Grid>(*rigidbody->get_grid()), rigidbody->get_waters(), rigidbody->fitter->fit_chi2_only());
}

Sequencer::~Sequencer() = default;

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

std::shared_ptr<fitter::Fit> Sequencer::execute() {
    for (auto& e : elements) {
        e->run();
    }
    rigidbody->update_fitter(rigidbody->fitter);
    return rigidbody->fitter->fit();
}