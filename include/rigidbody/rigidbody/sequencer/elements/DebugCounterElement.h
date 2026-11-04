// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class LoopElement;

    /**
     * @brief Test-only element which records each time it is run, so control flow can be checked without fitting.
     *
     * Script syntax: `debug_counter [name]`. It is a registered element but deliberately undocumented. A test finds the
     * instance in the parsed element tree after Sequencer::execute() and reads its state. The state lives on the
     * instance and is never reset, so a second execute() of the same Sequencer accumulates. It is not an optimization
     * step, so it does not change {iterations_total} or the global iteration counter.
     */
    class DebugCounterElement : public GenericElement {
        public:
            explicit DebugCounterElement(std::string name = "");
            ~DebugCounterElement() override;

            void run() override;

            /**
             * @brief The number of times this element has run.
             */
            int count() const;

            /**
             * @brief The global iteration counter (LoopElement::_get_current_iteration()) at each run, in order.
             */
            const std::vector<int>& iterations() const;

            const std::string& name() const;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            std::string label;
            std::vector<int> recorded;
    };
}
