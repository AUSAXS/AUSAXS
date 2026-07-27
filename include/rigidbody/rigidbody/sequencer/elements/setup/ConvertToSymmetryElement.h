// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <math/Vector3.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Collapse a set of loaded copy bodies into a single body carrying a fitted symmetry.
     *
     * The inverse of @ref SymmetryElement: fits the requested symmetry to the assembly formed by the given bodies (copies of the same molecule), 
     * installs it on the first (primary) body, and removes the redundant copies. The fit is rejected if its residual RMSD exceeds @ref tolerance.
     *
     * A setup-time operation. The number of bodies must equal repetitions()+1 for the symmetry (e.g. 4 for c4).
     *
     * As a special case, a single body is taken to be the assembly itself rather than one copy of it, and is split into the required number of copies before
     * the fit (see @ref _split_into_copies). This is the common case for a structure whose copies share one chain, which cannot be separated by chain id.
     */
    class ConvertToSymmetryElement : public GenericElement {
        public:
            //< Default residual-RMSD threshold (Å) above which the fit is considered a mismatch.
            static constexpr double default_tolerance = 2.0;

            /**
             * @param owner The owning sequencer.
             * @param bodies The participating body indices, primary first; the rest may be in any order. A single index requests the auto-split described above.
             * @param symmetry_name The target symmetry (e.g. "c4", "d3", "p2-p2").
             * @param tolerance Residual-RMSD threshold (Å); the setup fails if the fit exceeds it.
             */
            ConvertToSymmetryElement(observer_ptr<Sequencer> owner, std::vector<int> bodies, const std::string& symmetry_name, double tolerance = default_tolerance);
            ~ConvertToSymmetryElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            void _convert(const std::vector<int>& bodies, const std::string& symmetry_name, double tolerance);

            /**
             * @brief Decompose a single body that is itself the full assembly into the copies the symmetry expects.
             *
             * Splits the body's atoms into @p copies_wanted equal contiguous chunks and reduces the body itself (and its stored initial conformation) to the
             * leading chunk, leaving it as the reference copy the fitted symmetry regenerates the others from.
             *
             * @return The world-space coordinates of each copy, in file order, ready for the fit.
             * @throws sequencer::except::parse_error if the atom count does not divide evenly among the copies.
             */
            std::vector<std::vector<Vector3<double>>> _split_into_copies(int primary, std::size_t copies_wanted, const std::string& symmetry_name);

            observer_ptr<Sequencer> owner;
    };
}
