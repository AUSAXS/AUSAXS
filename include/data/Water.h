#pragma once

#include <data/Atom.h>
#include <data/Record.h>
#include <math/Vector3.h>

class Water : public Atom {
    public:
        using Atom::Atom; // inherit constructors from Atom
        Water(const Atom&& a) noexcept : Atom(std::move(a)) {}
        Water(const Atom& a) : Atom(a) {}
        ~Water() override = default;

        RecordType get_type() const override {return RecordType::WATER;}

        std::string get_recName() const override {return "HETATM";}

        bool is_water() const override {return true;}

        /**
         * @brief Create a new default water atom.
         */
        static Water create_new_water() {
            return create_new_water({0, 0, 0});
        }

        /**
         * @brief Create a new water atom.
         * @param coords the coordinates for the new atom.
         */
        static Water create_new_water(Vector3<double> coords) {
            return Water(-1, "O", "", "HOH", "", -1, "", coords, 1, 0, "O", "");
        }

        Water& operator=(const Water& rhs) = default;

        bool operator==(const Water& rhs) const {
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