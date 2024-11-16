#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <string>

namespace ausaxs::rigidbody::sequencer {
    class RelativeHydrationElement : public GenericElement {
        public:
            RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<double>& ratios);
            RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<double>& ratios, const std::vector<std::string>& names);
            ~RelativeHydrationElement() override;

            void run() override;

        private:
            observer_ptr<Sequencer> owner;
            std::vector<double> ratios;
    };
}