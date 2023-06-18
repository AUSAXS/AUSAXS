TEST_CASE("split_body", "[body],[files]") {
    settings::general::verbose = false;
    vector<int> splits = {9, 99};
    Protein protein = rigidbody::BodySplitter::split("test/files/LAR1-2.pdb", splits);

    // check sizes
    REQUIRE(protein.body_size() == 3);
    Body &b1 = protein.get_body(0), &b2 = protein.get_body(1), &b3 = protein.get_body(2);

    REQUIRE(b1.get_atoms().size() == 136);
    REQUIRE(b2.get_atoms().size() == 812-136);
    REQUIRE(b3.get_atoms().size() == 1606-812);

    // check start and end resseq
    CHECK(b1.get_atoms().back().resSeq == 8);
    CHECK(b2.get_atom(0).resSeq == 9);
    CHECK(b2.get_atoms().back().resSeq == 98);
    CHECK(b3.get_atom(0).resSeq == 99);
}
