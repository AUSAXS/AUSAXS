#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <data/DataFwd.h>

#include <rigidbody/sequencer/setup/SetupElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <utility/observer_ptr.h>

namespace rigidbody::sequencer {
    class Sequencer : public LoopElement, public SetupElement {
        public:
            Sequencer();
            Sequencer(const io::ExistingFile& saxs);
            ~Sequencer();

            std::shared_ptr<fitter::FitResult> execute() override;

            observer_ptr<RigidBody>& _get_rigidbody();
            observer_ptr<RigidBody> _get_rigidbody() const override;

            observer_ptr<detail::BestConf> _get_best_conf() const override;

            observer_ptr<const Sequencer> _get_sequencer() const override;

            bool _optimize_step() const;

        private:
            observer_ptr<RigidBody> rigidbody;
            std::unique_ptr<detail::BestConf> best;
    };
}