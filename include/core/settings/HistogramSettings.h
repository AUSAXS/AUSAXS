// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/ExvSettings.h>
#include <settings/SettingsHelper.h>

#include <string_view>

namespace ausaxs::settings {
    /// @brief Settings controlling the q-axis of the scattering curve and the distance-histogram binning.
    struct EXPORT axes {
        static int skip;                       // The number of points to skip from the top of the scattering curve.
        static detail::Setting<double> qmin;            // Lower limit on the used q-values
        static detail::Setting<double> qmax;            // Upper limit on the used q-values
        static detail::Setting<double> bin_width;       // The bin width to use for the distance histogram.
        static bool clamp_to_qrange;                    // Whether to clamp the input q-range to the range defined by qmin and qmax.
        static bool rebin;                              // Whether to rebin the scattering curve to increase the information content of each data point.
    };

    /// @brief Settings selecting how distance histograms are computed.
    struct EXPORT hist {
        /**
         * @brief The available kinds of histogram manager; see get_histogram_manager().
         *
         * A kind only decides how the histogram is recalculated. Which variant of it is used follows from the excluded volume
         * model, settings::exv::exv_method: the simple model weights each atom, the form factor models resolve the atoms by form
         * factor type, and the grid models have their own managers, which always recalculate everything.
         */
        enum class HistogramManagerChoice {
            HistogramManager,                    // A single-threaded manager that recalculates the entire histogram every time. A reference implementation, for the simple model only.
            HistogramManagerMT,                  // A multithreaded manager that recalculates the entire histogram every time.
            HistogramSymmetryManagerMT,          // A multithreaded manager for molecules with symmetries, evaluating each symmetric copy only once.
            PartialHistogramManager,             // A single-threaded manager that only recalculates the parts of the histogram that have been changed between each call. A reference implementation, for the simple model only.
            PartialHistogramManagerMT,           // A multithreaded implementation of the partial manager.
            PartialHistogramSymmetryManagerMT,   // A multithreaded implementation of the partial manager for molecules with symmetries.
            Count,
        };
        struct WeightedBins {
            enum class Value {
                True,
                False,
                Auto
            };

            WeightedBins() = default;
            WeightedBins(bool value) : value(value ? Value::True : Value::False) {}
            WeightedBins(std::string_view str);
            WeightedBins(WeightedBins::Value value) : value(value) {}

            bool is_auto() const {return value == Value::Auto;}
            bool is_true() const {return value == Value::True;}
            bool is_false() const {return value == Value::False;}

            Value value = Value::Auto;
        };
        static WeightedBins weighted_bins;

        /**
         * @brief Get the kind of histogram manager corresponding to the current number of threads and partial preference.
         */
        static HistogramManagerChoice get_histogram_manager();

        /**
         * @brief Check if a manager supports partial calculations, where only the contributions of a changed body are recalculated.
         *        These are the only managers suitable for iterative optimization, where a single body is moved between each evaluation.
         *
         * @param choice The kind of manager.
         * @param exv_method The excluded volume model, which decides the variant of the kind. The grid models have no partial variant.
         */
        static bool supports_partial_calculation(HistogramManagerChoice choice, exv::ExvMethod exv_method = exv::exv_method);
    };
}