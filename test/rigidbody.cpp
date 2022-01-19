#include "rigidbody/RigidTransform.h"
#include "rigidbody/RigidBody.h"
#include "data/Protein.h"

#include "catch2/catch.hpp"

// TEST_CASE("Constraints", "[rigidbody]") {
//     Atom a1(Vector3(-1, -1, -1), 1, "C", "C", 1);
//     Atom a2(Vector3(-1,  1, -1), 1, "C", "C", 1);
//     Atom a3(Vector3(-1, -1,  1), 1, "C", "C", 1);
//     Atom a4(Vector3(-1,  1,  1), 1, "C", "C", 1);
//     Atom a5(Vector3( 1, -1, -1), 1, "C", "C", 1);
//     Atom a6(Vector3( 1,  1, -1), 1, "C", "C", 1);
//     Atom a7(Vector3( 1, -1,  1), 1, "C", "C", 1);
//     Atom a8(Vector3( 1,  1,  1), 1, "He", "He", 1);

//     vector<Atom> b1 = {a1, a2};
//     vector<Atom> b2 = {a3, a4};
//     vector<Atom> b3 = {a5, a6};
//     vector<Atom> b4 = {a7, a8};
//     vector<vector<Atom>> ap = {b1, b2, b3, b4};
//     Protein protein(ap);

//     RigidBody rigidbody(protein);
    
//     SECTION("Invalid constraints") {
//         REQUIRE_THROWS(rigidbody.create_constraint(a1, a2)); // same body
//         REQUIRE_THROWS(rigidbody.create_constraint(a6, a8)); // non-C
//     }

//     SECTION("Check construction") {
//         rigidbody.create_constraint(a1, a3);

//         REQUIRE(*rigidbody.constraints[0].atom1 == a1);
//         REQUIRE(*rigidbody.constraints[0].atom2 == a3);
//         REQUIRE(rigidbody.constraints[0].body1->uid == protein.bodies[0].uid);
//         REQUIRE(rigidbody.constraints[0].body2->uid == protein.bodies[1].uid);
//     }
// }

TEST_CASE("RigidTransform", "[rigidbody]") {
    Atom a1(Vector3(-1, -1, -1), 1, "C", "C", 1);
    Atom a2(Vector3(-1,  1, -1), 1, "C", "C", 1);
    Atom a3(Vector3(-1, -1,  1), 1, "C", "C", 1);
    Atom a4(Vector3(-1,  1,  1), 1, "C", "C", 1);
    Atom a5(Vector3( 1, -1, -1), 1, "C", "C", 1);
    Atom a6(Vector3( 1,  1, -1), 1, "C", "C", 1);
    Atom a7(Vector3( 1, -1,  1), 1, "C", "C", 1);
    Atom a8(Vector3( 1,  1,  1), 1, "C", "C", 1);

    Body b1(std::vector<Atom>{a1, a2});
    Body b2(std::vector<Atom>{a3, a4});
    Body b3(std::vector<Atom>{a5, a6});
    Body b4(std::vector<Atom>{a7, a8});
    vector<Body> ap = {b1, b2, b3, b4};
    Protein protein(ap);

    RigidBody rigidbody(protein);

    Constraint constraint1(&a1, &a3, &b1, &b2);
    Constraint constraint2(&a5, &a7, &b3, &b4);
    rigidbody.add_constraint(constraint1);
    rigidbody.add_constraint(constraint2);

    RigidTransform transform(&rigidbody);

    // SECTION("get_connected") {
        vector<const Body*> group = transform.get_connected(constraint1);
        REQUIRE(group.size() == 1);
        REQUIRE(*group[0] == b1);

        group = transform.get_connected(constraint2);
        REQUIRE(group.size() == 1);
        REQUIRE(*group[0] == b3);

        Constraint constraint3(&a3, &a5, &b2, &b3);
        rigidbody.add_constraint(constraint3);
        group = transform.get_connected(constraint3);
        REQUIRE(group.size() == 2);
        REQUIRE(*group[0] == b1);
        REQUIRE(*group[1] == b2);

        Constraint constraint4(&a1, &a7, &b1, &b4);
        rigidbody.add_constraint(constraint4);
        group = transform.get_connected(constraint3);
        REQUIRE(group.size() == 4);
    // }
}