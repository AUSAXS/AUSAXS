// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/PhaseGuard.h>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/parse_error.h>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    std::string describe(const detail::ScriptPhase::Position& position) {
        return "\"" + position.element + "\"" + (position.line == 0 ? "" : " at line " + std::to_string(position.line));
    }

    // the error for an element appearing after a phase it must precede
    std::string too_late(const detail::ScriptPhase& state, Phase max) {
        // the first phase past the allowed ones that the script has entered; phases may be skipped, e.g. a loop directly after the load
        int closed = static_cast<int>(max) + 1;
        while (closed < static_cast<int>(state.current) && state.entered_by[closed].element.empty()) {++closed;}
        const auto& closer = state.entered_by[closed];

        switch (static_cast<Phase>(closed)) {
            case Phase::Structure:
                return "Only one \"load\" element is allowed per script, and the structure was already loaded by " + describe(closer) + ".";
            case Phase::Constraints:
                return "It must come before " + describe(closer) + ". The set of bodies cannot be changed once constraints have been "
                    "declared, as both the constraints and the symmetry targets are indexed by body.";
            case Phase::Run:
                return "It must come before " + describe(closer) + ", which starts the optimization. The setup cannot be changed afterwards.";
            default:
                throw ausaxs::except::unexpected("PhaseGuard: unknown phase.");
        }
    }
}

PhaseGuard::PhaseGuard(observer_ptr<LoopElement> owner, Phase min, Phase max, Phase advance) {
    if (owner == nullptr) {return;} // the root Sequencer is not part of the script
    auto& state = owner->_get_sequencer()->_script_phase();

    if (state.current < min) {
        throw except::parse_error(state.parsing.element, state.parsing.line, "It requires a loaded structure, so it must come after the \"load\" element.");
    }
    if (max < state.current) {
        throw except::parse_error(state.parsing.element, state.parsing.line, too_late(state, max));
    }
    if (state.current < advance) {
        state.current = advance;
        state.entered_by[static_cast<int>(advance)] = state.parsing;
    }
}
