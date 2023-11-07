#pragma once

#include <table/ArrayDebyeTable.h>
#include <table/DebyeTable.h>

namespace hist::detail {
    struct GenericDebyeTable {
        virtual 
    };

    struct DebyeTableSwitcher {
        template <typename T>
        static const table::DebyeTable& get_table() {
            return table::ArrayDebyeTable::get_instance();
        }
    };
}