// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody::sequencer {
    class GenericElement {
        public:
            GenericElement() = default;
            virtual ~GenericElement() = default;

            virtual void run() = 0;
    };
}