// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::rigidbody::sequencer {
    class ParameterElement : public GenericElement {
        public:
            ParameterElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy);
            ~ParameterElement() override;

            void run() override;

            observer_ptr<rigidbody::parameter::ParameterGenerationStrategy> get_parameter_strategy() const;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            observer_ptr<LoopElement> owner;
            std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy;
            int iterations = 0;
    };
}