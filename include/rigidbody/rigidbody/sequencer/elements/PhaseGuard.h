// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <array>
#include <string>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief The phases a script passes through, in the order they must appear in it.
     */
    enum class Phase {
        Empty,          // nothing is loaded yet
        Structure,      // the bodies are being established or reshaped
        Constraints,    // the set of bodies is fixed, and the constraints between them are being declared
        Run,            // the optimization has started
        Count
    };

    namespace detail {
        /**
         * @brief The phase a script has reached while it is being parsed, and which elements took it there.
         */
        struct ScriptPhase {
            struct Position {
                std::string element;
                int line = 0;
            };

            Phase current = Phase::Empty;
            Position parsing;                                                   // the element currently being parsed, stamped by the parser
            std::array<Position, static_cast<int>(Phase::Count)> entered_by;    // the element which first advanced the script to each phase
        };
    }

    /**
     * @brief Restricts where an element may appear in a script.
     *
     * Constructing a guard checks the phase the script has reached against the phases the element may appear in, and then
     * advances the script to the element's own phase. Since every element is constructed while its script is parsed, an
     * element out of place is a parse error, raised before the element does any work of its own.
     *
     * The named guards below are meant to be used in two ways: inherited by an element, so it can never be constructed out
     * of place, and constructed as a temporary at the top of an element's _parse function whenever that function itself needs
     * the structure. Checking twice is harmless for every guard which allows the phase it advances to, i.e. all but LoadPhaseGuard.
     */
    class PhaseGuard {
        protected:
            /**
             * @param owner The element owning the guarded one. The root Sequencer has none, and is never checked.
             * @param min The earliest phase the element may appear in.
             * @param max The latest phase the element may appear in.
             * @param advance The phase the element advances the script to. An element with a phase below the current one leaves it alone.
             */
            PhaseGuard(observer_ptr<LoopElement> owner, Phase min, Phase max, Phase advance);
    };

    /**
     * @brief Guard for the element loading the structure. Only a single one is allowed per script.
     */
    struct LoadPhaseGuard : PhaseGuard {
        explicit LoadPhaseGuard(observer_ptr<LoopElement> owner) : PhaseGuard(owner, Phase::Empty, Phase::Empty, Phase::Structure) {}
    };

    /**
     * @brief Guard for elements changing the set of bodies. These must appear before any constraints are declared, as both the
     *        constraints and the symmetry targets are indexed by body.
     */
    struct StructurePhaseGuard : PhaseGuard {
        explicit StructurePhaseGuard(observer_ptr<LoopElement> owner) : PhaseGuard(owner, Phase::Structure, Phase::Structure, Phase::Structure) {}
    };

    /**
     * @brief Guard for elements declaring constraints between the bodies. These must appear before the optimization starts.
     */
    struct ConstraintPhaseGuard : PhaseGuard {
        explicit ConstraintPhaseGuard(observer_ptr<LoopElement> owner) : PhaseGuard(owner, Phase::Structure, Phase::Constraints, Phase::Constraints) {}
    };

    /**
     * @brief Guard for elements which run the optimization. Once one appears, the setup is closed.
     */
    struct RunPhaseGuard : PhaseGuard {
        explicit RunPhaseGuard(observer_ptr<LoopElement> owner) : PhaseGuard(owner, Phase::Structure, Phase::Run, Phase::Run) {}
    };

    /**
     * @brief Guard for elements which only need a loaded structure. These may appear anywhere after it, and do not advance the phase.
     *        This is what lets e.g. a "parameters" element sit between the constraints and the loop it configures.
     */
    struct StructureRequiredGuard : PhaseGuard {
        explicit StructureRequiredGuard(observer_ptr<LoopElement> owner) : PhaseGuard(owner, Phase::Structure, Phase::Run, Phase::Empty) {}
    };
}
