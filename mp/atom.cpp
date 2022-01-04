#include "Test.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

void test_constructors() {
    // "element", "resName", and "name" are used for some internal logic, and must have reasonable values. "" can also be used. 
//*** ATOMS ***//
    Atom a1 = Atom(15, "", "altLoc", "", "chainID", 3, "iCode", Vector3({0, 1, 2}), 2.5, 3.5, "He", "2-");
    Atom a2 = Atom({3, 0, 5}, 2, "He", "", 3);

    IS_TRUE(a1.serial == 15);
    IS_TRUE(a1.name == "");
    IS_TRUE(a1.altLoc == "altLoc");
    IS_TRUE(a1.resName == "");
    IS_TRUE(a1.chainID == "chainID");
    IS_TRUE(a1.resSeq == 3);
    IS_TRUE(a1.iCode == "iCode");
    IS_TRUE(a1.coords == Vector3({0, 1, 2}));
    IS_TRUE(a1.occupancy == 2.5);
    IS_TRUE(a1.tempFactor == 3.5);
    IS_TRUE(a1.element == "He");
    IS_TRUE(a1.charge == "2-");

    IS_TRUE(a2.serial == 3);
    IS_TRUE(a2.name == "");
    IS_TRUE(a2.altLoc == "");
    IS_TRUE(a2.resName == "");
    IS_TRUE(a2.chainID == "");
    IS_TRUE(a2.resSeq == -1);
    IS_TRUE(a2.iCode == "");
    IS_TRUE(a2.coords == Vector3({3, 0, 5}));
    IS_TRUE(a2.occupancy == 2);
    IS_TRUE(a2.tempFactor == -1);
    IS_TRUE(a2.element == "He");
    IS_TRUE(a2.charge == "");

//*** HETATOMS ***//
    Hetatom w1 = Hetatom::create_new_water(Vector3({1, 2, 3}));
    IS_TRUE(w1.serial == -1);
    IS_TRUE(w1.name == "O");
    IS_TRUE(w1.altLoc == "");
    IS_TRUE(w1.resName == "HOH");
    IS_TRUE(w1.chainID == "");
    IS_TRUE(w1.resSeq == -1);
    IS_TRUE(w1.iCode == "");
    IS_TRUE(w1.coords == Vector3({1, 2, 3}));
    IS_TRUE(w1.occupancy == 1);
    IS_TRUE(w1.tempFactor == 0);
    IS_TRUE(w1.element == "O");
    IS_TRUE(w1.charge == "");
}

void test_operators() {
//*** ATOMS ***//
    Atom a1 = Atom({3, 0, 5}, 2, "He", "", 3);
    Atom a2 = a1;
    IS_TRUE(a1 == a2);

    a2 = Atom({0, 4, 1}, 2, "He", "", 2);
    IS_TRUE(a1 != a2);
    IS_TRUE(a2 < a1);

//*** HETATOMS ***//
    Hetatom w1 = Hetatom({3, 0, 5}, 2, "He", "", 3);
    Hetatom w2 = w1;
    IS_TRUE(w1 == w2);

    w2 = Atom({0, 4, 1}, 2, "He", "", 2);
    IS_TRUE(w1 != w2);
    IS_TRUE(w2 < w1);
}

int main(void) {
    cout << "Summary of Atom testing:" << endl;
    test_constructors();
    test_operators();
    if (passed_all) {
        cout << "\033[1;32m" << "All Atom tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Atom tests failed.");
    }
    return 0;
}