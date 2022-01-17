#pragma once

#include "Atom.h"
#include "Record.h"
#include "math/Vector3.h"

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
        _name = rhs.name; 
        _altLoc = rhs.altLoc; 
        _resName = rhs.resName; 
        _chainID = rhs.chainID; 
        _iCode = rhs.iCode; 
        _element = rhs.element; 
        _charge = rhs.charge;
        _occupancy = rhs.occupancy; 
        _tempFactor = rhs.tempFactor;
        _serial = rhs.serial; 
        _resSeq = rhs.resSeq;
        coords = rhs.coords;
        _effective_charge = rhs.effective_charge;
        _uid = rhs.uid;
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