// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef CONSTEXPR_LOOKUP_TABLE
    #include <table/ArrayDebyeTable.h>
    #include <settings/HistogramSettings.h>
    #include <utility/Logging.h>
    #include <utility/Utility.h>

    #include <string>

    using namespace ausaxs;
    using namespace ausaxs::table;

    void ArrayDebyeTable::check_default(const std::vector<double>& q, const std::vector<double>& d) {
        if (!logging::logging_enabled()) {return;}
        auto warn = [] (std::string_view reason) {
            logging::log("Warning in ArrayDebyeTable::check_default: Incompatible with default tables.\n\tReason: " + std::string(reason));
        };

        const Axis& axis = constants::axes::q_axis;
        auto qvals = axis.as_vector();
        int i = 0;
        for (; i < axis.bins; ++i) {
            if (utility::approx(q.front(), qvals[i])) {break;}
        }
        if (i == axis.bins) [[unlikely]] {
            warn("q[0] does not match any index of default q-array");
            return;
        }
        if (q[0] != qvals[i]) [[unlikely]] {
            warn("q[0] != axis.min");
        }

        if (q[1] != qvals[i+1]) [[unlikely]] {
            warn("q[1] != axis.min + (axis.max-axis.min)/axis.bins");
        }

        if (q[2] != qvals[i+2]) [[unlikely]] {
            warn("q[2] != axis.min + 2*(axis.max-axis.min)/axis.bins");
        }

        check_default(d);
    }

    void ArrayDebyeTable::check_default(const std::vector<double>& d) {
        if (!logging::logging_enabled()) {return;}
        auto warn = [] (std::string_view reason) {
            logging::log("Warning in ArrayDebyeTable::check_default: Incompatible with default tables.\n\tReason: " + std::string(reason));
        };

        // check empty
        if (d.empty()) [[unlikely]] {
            warn("d.empty()");
            return;
        }

        // check if too large for default table
        if (d.back() > constants::axes::d_axis.max) [[unlikely]] {
            warn("d.back() > default_size");
        }

        // check first width (d[1]-d[0] may be different from the default width)
        if (!utility::approx(d[2]-d[1], constants::axes::d_axis.width())) [[unlikely]] {
            warn("!utility::approx(d[2]-d[1], width)");
        }

        // check second width
        if (!utility::approx(d[3]-d[2], constants::axes::d_axis.width())) [[unlikely]] {
            warn("!utility::approx(d[3]-d[2], width)");
        }
    }

    #pragma message("Precompiling sinc lookup table. This may take a couple of minutes...")
    
    inline constexpr ArrayDebyeTable default_table;
    const ArrayDebyeTable& ArrayDebyeTable::get_default_table() {
        return default_table;
    }
#endif