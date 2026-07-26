// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/selection/ParameterMask.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <vector>

namespace ausaxs::rigidbody {
    namespace selection {
        /**
         * @brief This super-class defines the interface for the body selection strategies for the rigid-body optimization.
         * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
         */
        class BodySelectStrategy {
            public:
                BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody);
                virtual ~BodySelectStrategy() = default;

                /**
                 * @brief What a single optimization step acts on.
                 */
                struct Target {
                    unsigned int ibody = 0; //< the body to transform
                    int iconstraint = -1;   //< index into the body's constraint list, or -1 to transform the body on its own
                    int isymmetry = -1;     //< the body's symmetry slot to drive, or -1 to drive every one of them alike
                };

                struct SelectionResult {
                    unsigned int ibody;
                    int iconstraint;
                    int isymmetry;
                    ParameterMask mask;
                };

                /**
                 * @brief Get the next target to be transformed.
                 *
                 * @param mask The parameter mask this step will run under. Implementations must return a target whose move space is non-empty under it: in
                 *        particular, a symmetry-only mask requires a target naming a drivable symmetry slot (see symmetry_candidates).
                 */
                virtual Target next(const ParameterMask& mask) = 0;

                /**
                 * @brief Draw a parameter mask from the configured mask strategy, then a target compatible with it.
                 *
                 * Drawing in this order is what keeps the two halves consistent: the mask decides which parameter classes are live, and the target is then
                 * chosen so that at least one of them can actually move. Drawing them independently is how a step ends up perturbing nothing at all.
                 *
                 * @throws except::invalid_argument if the mask admits only symmetry parameters and the molecule declares no drivable symmetry.
                 */
                SelectionResult next_mask();

                /**
                 * @brief Replace the mask strategy used by next_mask(). Takes ownership.
                 */
                void set_mask_strategy(std::unique_ptr<ParameterMaskStrategy> strategy);

            protected:
                observer_ptr<const Rigidbody> rigidbody;

                /**
                 * @brief The current number of bodies in the molecule.
                 */
                unsigned int size_body() const;

                /**
                 * @brief True if the mask leaves the body's rigid pose frozen, so only a symmetry slot can move under it.
                 */
                static bool symmetry_only(const ParameterMask& mask);

                /**
                 * @brief Every symmetry slot the optimizer may drive. Non-owning views onto a shared symmetry are excluded, as driving them is a no-op.
                 */
                const std::vector<SymmetryTargets::Target>& symmetry_candidates() const;

                /**
                 * @brief Pick the constraint index for a body, mirroring the historical behaviour: none if it has no constraints, the only one if it has
                 *        exactly one, and a uniformly random one otherwise.
                 */
                int random_constraint(unsigned int ibody) const;

            private:
                std::unique_ptr<ParameterMaskStrategy> mask_strategy;
        };
    }
}
