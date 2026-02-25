// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/ParameterMask.h>
#include <utility/Random.h>

#include <random>

namespace ausaxs::rigidbody::selection {
    /**
     * @brief Abstract strategy that determines which transformation parameters are active each step.
     *
     * Implementations return a ParameterMask that is applied to the generated parameters
     * in BodySelectStrategy::next_masked() before they are passed to a transform strategy.
     */
    struct ParameterMaskStrategy {
        virtual ~ParameterMaskStrategy() = default;
        virtual ParameterMask next() = 0;
    };

    // Every parameter component is always active.
    struct AllMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override { return ParameterMask::all(); }
    };

    // Only the real body transform (translation + rotation) is active.
    struct RealOnlyMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override { return ParameterMask::real_only(); }
    };

    // Only symmetry parameters (offset, axis, angle) are active.
    struct SymmetryOnlyMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override { return ParameterMask::symmetry_only(); }
    };

    /**
     * @brief Alternates between real-transform step and symmetry-parameter step.
     *
     * Odd calls return real_only(); even calls return symmetry_only().
     * This ensures the optimizer dedicates equal attention to both parameter groups.
     */
    struct SequentialMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override {
            step_ = !step_;
            return step_ ? ParameterMask::real_only() : ParameterMask::symmetry_only();
        }
        private:
            bool step_ = false;
    };

    struct SequentialRealMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override {
            step_ = !step_;
            return step_ ? ParameterMask::real_only_rot() : ParameterMask::real_only_trans();
        }
        private:
            bool step_ = false;
    };

    struct SequentialSymmetryMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override {
            step_ = !step_;
            return step_ ? ParameterMask::symmetry_only() : ParameterMask::symmetry_only_axis();
        }
        private:
            int step_ = 0;
    };

    /**
     * @brief Picks either real-only or symmetry-only uniformly at random each step.
     */
    struct RandomMaskStrategy : ParameterMaskStrategy {
        ParameterMask next() override {
            static std::uniform_int_distribution<int> coin(0, 1);
            return coin(random::generator()) ? ParameterMask::real_only() : ParameterMask::symmetry_only();
        }
    };
}
