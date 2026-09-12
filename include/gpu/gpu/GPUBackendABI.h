// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstdint>
#include <type_traits>

/**
 * @brief The ABI between the library and a GPU backend.
 *
 * Backends are compiled separately (e.g. AdaptiveCpp for SYCL) and loaded via dlsym, so:
 * - everything has C linkage, only fixed-width trivially copyable types cross,
 * - nothing may throw (undefined across this boundary),
 * - nothing may be a template (a separate compile can't instantiate one),
 * - all memory is owned by the caller.
 */
namespace ausaxs::gpu::abi {
    extern "C" {
        /**
         * @brief The version of this interface.
         */
        constexpr std::int32_t version = 1;

        /**
         * @brief Status returned by the run functions.
         */
        enum Status : std::int32_t {
            ok = 0,
            no_device = 1,      // no usable device is present
            out_of_memory = 2,  // the device could not provide the buffers the jobs need
            kernel_failed = 3,  // the device rejected or failed to complete the work
            invalid_input = 4   // the jobs are malformed, e.g. a null coordinate pointer
        };

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
        static_assert(std::is_trivial_v<Job>, "Job must be trivial for ABI compatibility.");
        static_assert(std::is_standard_layout_v<Job>, "Job must be standard layout for ABI compatibility.");

        /**
         * @brief POD struct for a weighted histogram bin.
         */ 
        struct WeightedBin {
            double value;
            std::int64_t count;
            double center;
        };
        static_assert(std::is_trivial_v<WeightedBin>, "WeightedBin must be trivial for ABI compatibility.");
        static_assert(std::is_standard_layout_v<WeightedBin>, "WeightedBin must be standard layout for ABI compatibility.");

        /**
         * @brief The version of this interface the backend was compiled against.
         */
        using abi_version_fn = std::int32_t (*)();

        /**
         * @brief Whether a usable device is present. Never fails.
         */
        using available_fn = bool (*)();

        /**
         * @brief Name of the device the kernels run on, or "none" if there is no usable device.
         */
        using device_name_fn = const char* (*)();

        /**
         * @brief A message describing the most recent failure on the calling thread.
         *        Valid until the next call into the backend on that thread. Never null.
         */
        using last_error_fn = const char* (*)();

        /**
         * @brief Start a calculation, discarding anything left from a previous one.
         *
         * @param bin_count Number of bins per histogram. Distances beyond the last bin are discarded.
         * @param inv_width Inverse bin width; the bin of a distance d is round(inv_width*d).
         * @param weighted Whether the bins carry a count and a mean distance as well as a value.
         *
         * @return ok, or the reason the calculation could not be started. Never throws.
         */
        using begin_fn = Status (*)(std::int32_t bin_count, float inv_width, bool weighted);

        /**
         * @brief Queue correlations for evaluation. Returns without waiting for them.
         *        May be called any number of times between begin() and finish(). Jobs sharing a slot accumulate into the same histogram.
         *
         * @param jobs,n_jobs The correlations to evaluate. Their coordinates must stay alive until finish() returns.
         *
         * @return ok, or the reason the work could not be queued.
         */
        using submit_fn = Status (*)(const Job* jobs, std::int32_t n_jobs);

        /**
         * @brief Wait for everything submitted since begin() and read the histograms back.
         *        Only the slots that were actually submitted for hold a result. 
         *
         * @param n_slots Number of output histograms. 
         * @param out Caller-allocated output of n_slots*bin_count entries. Overwritten, not accumulated into.
         *
         * @return ok, or the reason the work could not be completed. Must match the @a weighted passed to begin(), or invalid_input is returned. Never throws.
         */
        using finish_unweighted_fn = Status (*)(std::int32_t n_slots, double* out);
        using finish_weighted_fn = Status (*)(std::int32_t n_slots, WeightedBin* out); //< @copydoc finish_unweighted_fn

        /**
         * @brief The entry points a backend must export.
         */
        std::int32_t ausaxs_gpu_abi_version();
        bool ausaxs_gpu_available();
        const char* ausaxs_gpu_device_name();
        const char* ausaxs_gpu_last_error();
        Status ausaxs_gpu_begin(std::int32_t bin_count, float inv_width, bool weighted);
        Status ausaxs_gpu_submit(const Job* jobs, std::int32_t n_jobs);
        Status ausaxs_gpu_finish_unweighted(std::int32_t n_slots, double* out);
        Status ausaxs_gpu_finish_weighted(std::int32_t n_slots, WeightedBin* out);

        // The symbol names the loader resolves, kept beside the declarations they must match.
        constexpr const char* symbol_abi_version       = "ausaxs_gpu_abi_version";
        constexpr const char* symbol_available         = "ausaxs_gpu_available";
        constexpr const char* symbol_device_name       = "ausaxs_gpu_device_name";
        constexpr const char* symbol_last_error        = "ausaxs_gpu_last_error";
        constexpr const char* symbol_begin             = "ausaxs_gpu_begin";
        constexpr const char* symbol_submit            = "ausaxs_gpu_submit";
        constexpr const char* symbol_finish_unweighted = "ausaxs_gpu_finish_unweighted";
        constexpr const char* symbol_finish_weighted   = "ausaxs_gpu_finish_weighted";
    }

    // sanity checks that the backend's symbols match the expected types
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_abi_version), abi_version_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_available), available_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_device_name), device_name_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_last_error), last_error_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_begin), begin_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_submit), submit_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_finish_unweighted), finish_unweighted_fn>);
    static_assert(std::is_same_v<decltype(&ausaxs_gpu_finish_weighted), finish_weighted_fn>);
}
