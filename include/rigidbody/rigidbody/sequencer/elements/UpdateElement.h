// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Publishes the current structure to a shared buffer during a run so the GUI can poll it and live-update its preview. 
     */
    class UpdateElement : public LoopElementCallback, public GenericElement {
        public:
            UpdateElement(observer_ptr<LoopElement> owner);
            ~UpdateElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

            static void lock();
            static void unlock();

            // latest published structure (explicit, symmetries realized), guarded by `mutex`. version is bumped on each publish and reset to 0 when a new sequence is parsed.
            inline static std::vector<double> x, y, z;
            inline static int version = 0;

            // set true by a consumer (e.g. the GUI) that polls the live structure. If false, the element publishes nothing.
            inline static bool live_consumer_connected = false;
    };
}
