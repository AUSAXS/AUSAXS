#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>

#include <algorithm>
#include <random>
#include <vector>

namespace test::form_factor {
    /**
     * @brief Activate a random form factor selection of (at most) n types, and return it.
     *        As in production, EXCLUDED_VOLUME and WATER come first and OTHER last; the n-3 types in between are drawn at random and in random order,
     *        so a test cannot accidentally rely on the active slot of a type matching its enum index.
     *
     *        This is only intended for tests of the product tables themselves.
     *        Tests computing a histogram of a molecule should select for that molecule with form_factor::manager::use_form_factors, like the histogram factory does.
     */
    inline std::vector<int> use_random_form_factors(int n = 10) {
        using ausaxs::form_factor::form_factor_t;
        constexpr int exv   = static_cast<int>(form_factor_t::EXCLUDED_VOLUME);
        constexpr int water = static_cast<int>(form_factor_t::WATER);
        constexpr int other = static_cast<int>(form_factor_t::OTHER);

        std::vector<int> pool;
        for (int i = 0; i < ausaxs::form_factor::total_ff_count; ++i) {
            if (i == exv || i == water || i == other) {continue;}
            pool.push_back(i);
        }
        std::mt19937 rng(std::random_device{}());
        std::ranges::shuffle(pool, rng);
        pool.resize(std::clamp<int>(n-3, 0, static_cast<int>(pool.size())));

        std::vector<int> selection{exv, water};
        selection.insert(selection.end(), pool.begin(), pool.end());
        selection.push_back(other);
        ausaxs::form_factor::manager::detail::use_form_factors(selection);
        return selection;
    }
}
