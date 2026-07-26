// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief Reject an element that would change the set of bodies once constraints have been declared.
     *
     * Constraints and the symmetry target pool are both indexed by body index and are rebuilt only while the setup phase is still running; adding or removing
     * bodies afterwards silently invalidates both. Declaring a constraint also requires naming the bodies it joins, so any script that mutates the structure
     * afterwards is describing a body it has already tied down - user error rather than a case worth supporting.
     *
     * @param element The element name to attribute the error to.
     * @throws sequencer::except::parse_error if any discoverable constraint has been declared.
     */
    void require_mutable_structure(observer_ptr<Sequencer> owner, const std::string& element);

    /**
     * @brief Erase the bodies at the given indices from the molecule and its conformation, and update the
     *        setup's body-name map so surviving names are reindexed to their shifted positions while the
     *        names of the erased bodies are dropped.
     *
     * @param indices Indices of the bodies to erase. May be given in any order; must not contain duplicates.
     */
    void erase_bodies(observer_ptr<Sequencer> owner, std::vector<int> indices);
}
