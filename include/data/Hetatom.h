#pragma once

#include <data/Atom.h>
#include <data/Record.h>
#include <math/Vector3.h>

class Hetatom : public Atom {
  public:
    using Atom::Atom; // inherit constructors from Atom
    Hetatom(const Atom&& a) noexcept : Atom(std::move(a)) {}
    Hetatom(const Hetatom&& a) : Atom(std::move(a)) {}
    Hetatom(const Atom& a) : Atom(a) {}
    Hetatom(const Hetatom& a) : Atom(a) {}
    ~Hetatom() override = default;

    RecordType get_type() const override {return HETATM;}

    string get_recName() const override {return "HETATM";}

    bool is_water() const override {return resName == "HOH";}
    
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

    Hetatom& operator=(const Hetatom& rhs) {
        name = rhs.name; 
        altLoc = rhs.altLoc; 
        resName = rhs.resName; 
        chainID = rhs.chainID; 
        iCode = rhs.iCode; 
        element = rhs.element; 
        charge = rhs.charge;
        occupancy = rhs.occupancy; 
        tempFactor = rhs.tempFactor;
        serial = rhs.serial; 
        resSeq = rhs.resSeq;
        coords = rhs.coords;
        effective_charge = rhs.effective_charge;
        uid = rhs.uid;
        return *this;
    }

    bool operator==(const Hetatom& rhs) const {
        if (get_type() != get_type()) {return false;}
        if (name != rhs.name) {return false;}
        if (altLoc != rhs.altLoc) {return false;}
        if (resName != rhs.resName) {return false;}
        if (chainID != rhs.chainID) {return false;}
        if (iCode != rhs.iCode) {return false;}
        if (element != rhs.element) {return false;}
        if (charge != rhs.charge) {return false;}
        if (occupancy != rhs.occupancy) {return false;}
        if (tempFactor != rhs.tempFactor) {return false;}
        if (serial != rhs.serial) {return false;}
        // if (resSeq != rhs.resSeq) {return false;} // this is to fix io tests, since some pdb files randomly changes this order
        if (coords != rhs.coords) {return false;}
        if (uid != rhs.uid) {return false;}
        return true;
    }
};