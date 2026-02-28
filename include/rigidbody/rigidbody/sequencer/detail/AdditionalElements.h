// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/Sequencer.h>

namespace ausaxs::rigidbody::sequencer::detail {
    struct SeedElement {static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);};
    struct LoopEndElement {static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);};
    struct OverlapStrengthElement {static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);};
    struct LogElement {static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);};
}