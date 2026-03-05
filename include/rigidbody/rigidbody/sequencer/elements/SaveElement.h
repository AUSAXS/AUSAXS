// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <utility/observer_ptr.h>
#include <io/File.h>

namespace ausaxs::rigidbody::sequencer {
    class SaveElement : public LoopElementCallback, public GenericElement {
        public:
            SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path);
            ~SaveElement() override;

            void run() override;

            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            io::File path;
    };
}