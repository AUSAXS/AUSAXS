// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Split an existing body into several new bodies at the given residue sequence ids.
     *
     * If the body carries one or more symmetries, each is turned into a ReferenceSymmetry shared by every resulting fragment with the first fragment becomes
     * the owning primary, and the rest hold non-owning views.
     */
    class SplitElement : public GenericElement {
        public:
            /**
             * @param owner       The owning sequencer.
             * @param body_name   Name of the (base) body to split.
             * @param splits      Residue sequence ids to split at; produces splits.size()+1 fragments.
             */
            SplitElement(observer_ptr<Sequencer> owner, const std::string& body_name, std::vector<int> splits);
            ~SplitElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            void _split(const std::string& body_name, const std::vector<int>& splits);
            observer_ptr<Sequencer> owner;
    };
}
