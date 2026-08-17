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

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

            /**
             * @brief Discard the output state shared by every save element: the '%' file counters and the open
             *        trajectory writers.
             *
             * That state is deliberately process-wide so a single script numbers its files consistently and appends to
             * one trajectory throughout. A host that runs several refinements in one process must call this between
             * them, or the next run continues the previous run's numbering and keeps writing to its trajectory file —
             * which, if that file has since been moved, means writing into the old location.
             */
            static void _reset_output_state();

        private:
            io::File path;
    };
}