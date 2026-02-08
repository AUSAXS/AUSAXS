#include <catch2/catch_test_macros.hpp>

#include <rigidbody/transform/BackupBody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::transform;

TEST_CASE("BackupBody::construction") {
    data::Molecule molecule("tests/files/2epe.pdb");
    const auto& body = molecule.get_body(0);
    
    SECTION("construct with body and index") {
        BackupBody backup(body, 0);
        
        CHECK(backup.index == 0);
        CHECK(backup.body.size_atom() == body.size_atom());
    }

    SECTION("construct with different index") {
        BackupBody backup(body, 5);
        
        CHECK(backup.index == 5);
    }
}

TEST_CASE("BackupBody::basic usage") {
    data::Molecule molecule("tests/files/2epe.pdb");
    auto& body = molecule.get_body(0);
    
    SECTION("stores independent copy of body") {
        BackupBody backup(body, 0);
        auto original_pos = body.get_atom(0).coordinates();
        
        body.translate(Vector3<double>(10, 10, 10));
        
        CHECK(backup.body.get_atom(0).coordinates() == original_pos);
        CHECK(body.get_atom(0).coordinates() != original_pos);
    }

    SECTION("can restore from backup") {
        BackupBody backup(body, 0);
        auto original_pos = body.get_atom(0).coordinates();
        
        body.translate(Vector3<double>(10, 10, 10));
        body = backup.body;
        
        CHECK(body.get_atom(0).coordinates() == original_pos);
    }
}
