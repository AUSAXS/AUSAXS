// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/WebGPU.h>

namespace ausaxs::gpu {
    /**
     * @brief The WebGPU instance, adapter and device used for all calculations.
     *
     * Acquiring a device is expensive (tens to hundreds of milliseconds), so a single instance is
     * shared by all backends through get(). It is created on first use and lives until the process
     * exits; the handles are deliberately never released, since their teardown order relative to the
     * backend implementation is not guaranteed.
     */
    class GPUInstance {
        public:
            /**
             * @brief Get the shared instance, creating it on first use.
             * @throws except::runtime_error if no device could be acquired.
             */
            static GPUInstance& get();

            /**
             * @brief Check whether a device can be acquired on this machine.
             *        Unlike get(), this does not throw.
             */
            static bool available();

            wgpu::Instance instance;
            wgpu::Device device;
            wgpu::Queue queue;

            /// @brief Let the WebGPU implementation run its pending callbacks.
            void process();

            /// @brief Block until @a done is set by a callback.
            void wait(bool& done);

        private:
            GPUInstance();
    };
}
