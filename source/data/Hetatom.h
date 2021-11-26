#pragma once

#include "Atom.h"
#include "Record.h"
#include "math/Vector3.h"

class Hetatom : public Atom {
public:
    using Atom::Atom; // inherit constructors from Atom
    ~Hetatom() override {}

    RecordType get_type() const override {return HETATM;}

    string get_recName() const override {return "HETATM";}

    bool is_water() const override {
        if (resName == "HOH") {
            return true;
        }
        return false;
    }

    /**
     * @brief Create a new water Atom.
     * @return A pointer to the new water Atom. 
     */
    static shared_ptr<Hetatom> create_new_water() {
        return create_new_water({0, 0, 0});
    }

    /**
     * @brief Create a new water Atom.
     * @param coords the coordinates for the new Atom.
     * @return A pointer to the new water Atom. 
     */
    static shared_ptr<Hetatom> create_new_water(Vector3 coords) {
        return std::make_shared<Hetatom>(-1, "O", "", "HOH", "", -1, "", coords, 1, 0, "O", "");
    }
};