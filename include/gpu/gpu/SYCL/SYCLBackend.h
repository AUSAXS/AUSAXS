// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstdint>

/**
 * @brief Narrow interface to the SYCL distance histogram kernels.
 *
 * The implementation is compiled by the SYCL compiler (AdaptiveCpp) rather than the compiler used
 * for the rest of the project, so nothing SYCL-specific may appear here. Only trivially copyable
 * types cross this boundary, and all memory is owned by the caller.
 */
namespace ausaxs::gpu::sycl_backend {
    /**
     * @brief The scale of the fixed-point bin values.
     *
     * The kernels cannot accumulate in floating point without losing most of the contributions to
     * rounding — the bins reach ~1e11, where the f32 spacing is larger than a typical contribution.
     * Bins are therefore 64-bit fixed-point integers, which is also exact and order-independent, so
     * the results are reproducible between runs. This matches the WebGPU backend, which needs the
     * same treatment for the same reason but has to emulate the 64-bit adds in two 32-bit words.
     */
    constexpr double fixed_point_scale = 65536;

    /**
     * @brief One queued correlation.
     *
     * The coordinates are read as tightly packed [x, y, z, w] floats, i.e. the layout of
     * hist::detail::CompactCoordinates. They must stay alive until the run call returns.
     */
    struct Job {
        const float* a1;        // first coordinate set
        const float* a2;        // second coordinate set, or nullptr for a self-correlation
        std::uint32_t n1, n2;   // number of atoms in each set; n2 is ignored for self-correlations
        std::uint32_t scaling;  // multiplicative factor applied to every contribution
        std::uint32_t slot;     // index of the output histogram this job accumulates into
    };

    /// @brief A weighted bin as written by the kernels. Value and center are fixed-point, count is exact.
    struct WeightedBin {
        std::int64_t value;
        std::int64_t count;
        std::int64_t center;
    };

    /**
     * @brief Check whether a usable device is present. Never throws.
     *        Calling this also creates the shared queue, so the cost of doing so is not attributed
     *        to the first calculation.
     */
    bool available();

    /// @brief Name of the device the kernels run on, or "none" if there is no usable device.
    const char* device_name();

    /**
     * @brief Evaluate all queued correlations.
     *
     * @param jobs,n_jobs The correlations to evaluate. Jobs sharing a slot accumulate into the same histogram.
     * @param inv_width Inverse bin width; the bin of a distance d is round(inv_width*d).
     * @param bin_count Number of bins per histogram. Distances beyond the last bin are discarded.
     * @param n_slots Number of output histograms.
     * @param out Caller-allocated output of n_slots*bin_count entries. Overwritten, not accumulated into.
     *
     * @throws std::runtime_error if no device is available or a kernel fails.
     */
    void run_unweighted(
        const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, std::int64_t* out
    );
    void run_weighted( //< @copydoc run_unweighted
        const Job* jobs, int n_jobs, float inv_width, int bin_count, int n_slots, WeightedBin* out
    );
}
