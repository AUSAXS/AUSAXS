#include "rigidbody/RigidTransform.h"
#include "rigidbody/RigidBody.h"
#include "data/Protein.h"

#include "catch2/catch.hpp"

TEST_CASE("Constraints", "[rigidbody]") {
    std::shared_ptr<Atom> a1 = std::make_shared<Atom>(Vector3(-1, -1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a2 = std::make_shared<Atom>(Vector3(-1,  1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a3 = std::make_shared<Atom>(Vector3(-1, -1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a4 = std::make_shared<Atom>(Vector3(-1,  1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a5 = std::make_shared<Atom>(Vector3( 1, -1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a6 = std::make_shared<Atom>(Vector3( 1,  1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a7 = std::make_shared<Atom>(Vector3( 1, -1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a8 = std::make_shared<Atom>(Vector3( 1,  1,  1), 1, "C", "C", 1);

    vector<Atom> b1 = {*a1, *a2};
    vector<Atom> b2 = {*a3, *a4};
    vector<Atom> b3 = {*a5, *a6};
    vector<Atom> b4 = {*a7, *a8};
    vector<vector<Atom>> ap = {b1, b2, b3, b4};
    Protein protein(ap, {});

    RigidBody rigidbody(protein);
    rigidbody.create_constraint(a1, a2);
    
}

TEST_CASE("RigidTransform", "[rigidbody]") {
    std::shared_ptr<Atom> a1 = std::make_shared<Atom>(Vector3(-1, -1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a2 = std::make_shared<Atom>(Vector3(-1,  1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a3 = std::make_shared<Atom>(Vector3(-1, -1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a4 = std::make_shared<Atom>(Vector3(-1,  1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a5 = std::make_shared<Atom>(Vector3( 1, -1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a6 = std::make_shared<Atom>(Vector3( 1,  1, -1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a7 = std::make_shared<Atom>(Vector3( 1, -1,  1), 1, "C", "C", 1);
    std::shared_ptr<Atom> a8 = std::make_shared<Atom>(Vector3( 1,  1,  1), 1, "C", "C", 1);

    vector<Atom> b1 = {*a1, *a2};
    vector<Atom> b2 = {*a3, *a4};
    vector<Atom> b3 = {*a5, *a6};
    vector<Atom> b4 = {*a7, *a8};
    vector<vector<Atom>> ap = {b1, b2, b3, b4};
    Protein protein(ap, {});

    RigidBody rigidbody(protein);
    rigidbody.create_constraint(a1, a2);

    RigidTransform transform(&rigidbody);
}