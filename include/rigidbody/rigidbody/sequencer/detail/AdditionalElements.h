// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/Sequencer.h>

namespace ausaxs::rigidbody::sequencer::detail {
    struct SeedElement {
        static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
        static std::vector<std::string> _valid_arguments();
    };
    struct LoopEndElement {
        static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
        static std::vector<std::string> _valid_arguments();
    };
    struct OverlapStrengthElement {
        static void _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
        static std::vector<std::string> _valid_arguments();
    };
    struct LogElement {
        static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
        static std::vector<std::string> _valid_arguments();
    };
}