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

            /**
             * @brief Define a single reference symmetry shared across several bodies.
             *
             * The first body owns the (cyclic) symmetry; the remaining bodies are linked to it
             * through non-owning views, so the whole group is replicated as one rigid assembly.
             */
            SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& body_names, const std::string& reference_symmetry);
            ~SymmetryElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            /**
             * @brief Install one pre-built symmetry per body, wiring up names, parameters and the grid.
             */
            void _add(const std::vector<std::string>& names, std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetries);

            /**
             * @brief Install a single reference symmetry shared across the named bodies.
             */
            void _add_reference(const std::vector<std::string>& body_names, const std::string& reference_symmetry);

            observer_ptr<Sequencer> owner;
    };
}