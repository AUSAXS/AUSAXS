// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <utility/observer_ptr.h>

#include <mutex>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Publishes the current structure to a shared buffer during a run so the GUI can poll
     *        it and live-update its preview. Only `update structure` is supported for now.
     *
     * run() fires on the optimisation (worker) thread while the GUI polls from another thread, so
     * the buffer is mutex-guarded and carries a version counter the reader uses to detect changes.
     * Only coordinates are published: the atom ordering is fixed after the setup phase, so the GUI
     * maps the backbone once (via preview_structure) and reuses that mask for every live frame.
     */
    class UpdateElement : public LoopElementCallback, public GenericElement {
        public:
            UpdateElement(observer_ptr<LoopElement> owner);
            ~UpdateElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

            // latest published structure (explicit, symmetries realized), guarded by `mutex`.
            // version is bumped on each publish and reset to 0 when a new sequence is parsed.
            inline static std::mutex mutex;
            inline static std::vector<double> x, y, z;
            inline static int version = 0;
    };
}
