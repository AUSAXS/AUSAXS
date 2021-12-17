#pragma once

#include "Atom.h"
#include "Record.h"
#include "math/Vector3.h"

class Hetatom : public Atom {
public:
    using Atom::Atom; // inherit constructors from Atom
    Hetatom(const Atom& a) : Atom(std::move(a)) {}
    // Hetatom(const Hetatom& a) : Atom(std::move(a)) {}
    ~Hetatom() override {}

    RecordType get_type() const override {return HETATM;}

    string get_recName() const override {return "HETATM";}

    bool is_water() const override {
        cout << "CHECKING HETATOM: " << resName << endl;
        return resName == "HOH";
    }
    
    /**
     * @brief Create a new default water atom.
     */
    static Hetatom create_new_water() {
        return create_new_water({0, 0, 0});
    }

    /**
     * @brief Create a new water atom.
     * @param coords the coordinates for the new atom.
     */
    static Hetatom create_new_water(Vector3 coords) {
        return Hetatom(-1, "O", "", "HOH", "", -1, "", coords, 1, 0, "O", "");
    }
};