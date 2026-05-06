#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorManager::set_custom_form_factors") {
    data::Molecule molecule("tests/files/2epe.pdb");
    FormFactorManager::set_custom_form_factors(molecule);
    volatile auto custom_tables = FormFactorManager::get_custom_tables();
}