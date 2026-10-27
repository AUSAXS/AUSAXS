// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef CONSTEXPR_LOOKUP_TABLE
    #include <table/ArrayDebyeTable.h>

    using namespace ausaxs;
    using namespace ausaxs::table;

    #pragma message("Precompiling sinc lookup table. This may take a couple of minutes...")
    
    inline constexpr ArrayDebyeTable default_table;
    const ArrayDebyeTable& ArrayDebyeTable::get_default_table() {
        return default_table;
    }
#endif