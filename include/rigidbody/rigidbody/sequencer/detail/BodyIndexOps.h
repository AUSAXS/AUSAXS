// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief Erase the bodies at the given indices from the molecule and its conformation, and update the
     *        setup's body-name map so surviving names are reindexed to their shifted positions while the
     *        names of the erased bodies are dropped.
     *
     * @param indices Indices of the bodies to erase. May be given in any order; must not contain duplicates.
     */
    void erase_bodies(observer_ptr<Sequencer> owner, std::vector<int> indices);
}
