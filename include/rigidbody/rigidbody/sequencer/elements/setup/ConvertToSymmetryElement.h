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
     * @brief Collapse a set of already-loaded copy bodies into a single reference body carrying a
     *        fitted symmetry.
     *
     * This is the inverse of expanding one body + a symmetry into many (@ref ConvertSymmetryElement):
     * given N bodies that are copies of the same molecule assembled with some point-group symmetry,
     * it fits the parameters of the requested symmetry to that assembly, installs the fitted symmetry
     * on the first (primary) body, and removes the now-redundant copies. The residual RMSD of the fit
     * is reported so the user can judge how well the symmetry actually holds.
     *
     * This is a setup-time operation. The number of participating bodies must equal
     * repetitions()+1 for the requested symmetry (e.g. exactly 4 bodies for c4); otherwise the setup
     * fails. Bodies must be given in symmetry (repetition) order.
     */
    class ConvertToSymmetryElement : public GenericElement {
        public:
            /**
             * @param owner The owning sequencer.
             * @param bodies The participating body indices, primary first, in repetition order.
             * @param symmetry_name The target symmetry (e.g. "c4", "d3", "tetrahedral").
             */
            ConvertToSymmetryElement(observer_ptr<Sequencer> owner, std::vector<int> bodies, const std::string& symmetry_name);
            ~ConvertToSymmetryElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            void _convert(const std::vector<int>& bodies, const std::string& symmetry_name);
            observer_ptr<Sequencer> owner;
    };
}
