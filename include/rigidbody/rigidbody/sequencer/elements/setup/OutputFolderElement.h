// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>
#include <io/IOFwd.h>

namespace ausaxs::rigidbody::sequencer {
    class OutputFolderElement : public GenericElement {
        public:
            /**
             * @brief Different modes for defining the output folder.
             */
            enum class Mode {
                RELATIVE_TERMINAL,
                RELATIVE_CONFIG,
            };

            OutputFolderElement(observer_ptr<Sequencer> owner, const io::Folder& path, Mode mode = Mode::RELATIVE_TERMINAL);
            ~OutputFolderElement() override = default;

            void run() override;

        private: 
            observer_ptr<Sequencer> owner;
    };
}