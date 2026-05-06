#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorManager::set_custom_form_factors") {
    data::Molecule molecule("tests/files/2epe.pdb");
    std::vector<data::AtomFF> atoms(molecule.get_body(0).get_atoms().begin(), molecule.get_body(0).get_atoms().begin()+10);
    int counter = 0;
    for (auto& a : atoms) {
        std::cout << "atoms[" << counter++ << "] = " << form_factor::to_string(a.form_factor_type()) << std::endl;
    }
    counter = 0;
    FormFactorManager::set_custom_form_factors(molecule);
    for (auto& a : atoms) {
        std::cout << "atoms[" << counter++ << "] = " << form_factor::to_string(a.form_factor_type()) << std::endl;
    }

    // auto custom_tables = FormFactorManager::get_custom_tables();
}