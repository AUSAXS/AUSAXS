// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody::sequencer {
    class SequenceElement {
        public:
            virtual void execute() = 0;
    };
}