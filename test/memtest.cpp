#include "catch2/catch.hpp"

#include "data/Atom.h"
#include "data/Protein.h"
#include "data/BodySplitter.h"

/**
 * @brief These tests are meant to be run with valgrind to help identify memory issues. 
 *        I'll expand it as I encounter more leaks/bugs with the tool. 
 */
TEST_CASE("atom_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);
    b2 = b1;
}

TEST_CASE("body_splitter", "[memtest]") {
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2.pdb", splits);
    RigidBody body(protein);
}

TEST_CASE("body_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);

    b2 = b1;
}

TEST_CASE("grid_add", "[memtest]") {
    vector<int> splits = {9, 99};
    Protein protein(BodySplitter::split("data/LAR1-2.pdb", splits));
    RigidBody rbody(protein);
    rbody.protein.generate_new_hydration();
    rbody.protein.clear_grid();
    rbody.protein.generate_new_hydration();
}

TEST_CASE("compact_coordinates", "[memtest]") {
    
}

TEST_CASE("file_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    vector<Hetatom> waters;

    std::shared_ptr<File> file2 = std::make_shared<File>();
    std::shared_ptr<File> file1 = std::make_shared<File>(atoms, waters);
    vector<Atom>& file2a = file2->protein_atoms;
    file2 = std::make_shared<File>(*file1);
    file2a = file2->protein_atoms;

    std::cout << file2a[0].as_pdb() << std::endl;
}

// TEST_CASE("body_splitter", "[memtest]") {
//     vector<int> splits = {9, 99};

//     Body body("data/LAR1-2.pdb");
//     vector<Atom>& atoms = body.protein_atoms;

//     // we define a boolean vector with one entry for each residue sequence id
//     vector<bool> split(atoms.back().resSeq, false);

//     // we then mark the ids where we want to split as true
//     std::for_each(splits.begin(), splits.end(), [&split] (const int id) {split[id] = true;});

//     vector<Body> bodies(splits.size()+1);
//     int index_body = 0; // current index in the bodies vector

//     // the two iterators marks the indices in the atoms vector where we want to split next time
//     vector<Atom>::const_iterator begin = atoms.begin(); // start at the beginning
//     vector<Atom>::const_iterator end; // no defined end yet
//     for (unsigned int i = 0; i < atoms.size(); i++) {
//         int resSeq = std::max(atoms[i].resSeq, 0); // in some files resSeq starts negative

//         // we can now in constant time look in our split vector to see if we should split at this atom
//         if (split[resSeq]) {
//             end = atoms.begin() + i;        // define the end index
//             vector<Atom> a(begin, end);     // create a new vector of atoms based on the start and end iterators
//             bodies[index_body++] = Body(a); // create a body from this vector
//             begin = end;                    // change the start index for the next split
//             split[resSeq] = false;          // mark it as false so we won't split again on the next atom
//         }
//     }

//     // // add the final body
//     // begin = end;
//     // end = atoms.end();
//     // vector<Atom> a(begin, end);
//     // bodies[index_body] = Body(a);

//     // Protein protein(bodies);
// }