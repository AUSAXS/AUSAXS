// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Define a symmetry for a given body.
     */
    class SymmetryElement : public GenericElement {
        public:
            SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry);

            /**
             * @brief Define symmetries given by name string, one per body.
             *
             * Names may be composite (e.g. "p2-c3"); see symmetry::create.
             */
            SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<std::string>& symmetry_names);
            ~SymmetryElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            /**
             * @brief Install one pre-built symmetry per body, wiring up names, parameters and the grid.
             */
            void _add(const std::vector<std::string>& names, std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetries);

            observer_ptr<Sequencer> owner;
    };
}