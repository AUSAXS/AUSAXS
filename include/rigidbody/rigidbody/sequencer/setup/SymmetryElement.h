#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Define a symmetry for a given body.
     */
    class SymmetryElement : public GenericElement {
        public:
            SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry);
            ~SymmetryElement() override;

            void run() override;

        private:
            observer_ptr<Sequencer> owner;
    };
}