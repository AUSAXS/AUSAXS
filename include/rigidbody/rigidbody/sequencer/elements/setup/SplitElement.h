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
     * The in-memory analogue of the load-time "split" (see BodySplitter): rather than splitting a freshly-read
     * PDB file, this partitions a body already present in the molecule (e.g. one produced by convert_to_symmetry)
     * using its per-atom residue-sequence metadata, which requires settings::molecule::store_residue_seq - already
     * enabled unconditionally by the sequencer.
     *
     * If the body carries one or more symmetries, each is turned into a ReferenceSymmetry shared by every
     * resulting fragment (the first fragment becomes the owning primary, the rest hold non-owning views), so the
     * fragments stay tied to the same symmetric assembly instead of each getting an independently-optimizable copy.
     *
     * A setup-time operation; not supported on a body that already participates in a symmetry shared with bodies
     * outside the split (i.e. a symmetry replica, or a body/view involved in another body's ReferenceSymmetry).
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
