#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorManager::set_custom_form_factors") {
    data::Molecule molecule("tests/files/2epe.pdb");
    std::vector<form_factor_t> original_ff_types, new_ff_types;

    int counter = 1;
    for (auto& a : molecule.iterate_atoms()) {
        original_ff_types.emplace_back(a.form_factor_type());
        if (10 < ++counter) {break;}
    }

    FormFactorManager::set_custom_form_factors(molecule);
}