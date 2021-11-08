#pragma once

// includes
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <boost/format.hpp>
#include <cstdio>

// ROOT
#include <TVector3.h>

// my own stuff
#include "Record.h"
#include "Tools.cpp"

using namespace ROOT;
using std::vector, std::string, std::cout, std::endl, std::setw;
using boost::format;

class Atom : public Record {
public:
    /** 
     * @brief Constructor for the Atom class. 
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom
     * @param symbol the atomic symbol of the base atom
     * @param comp the molecule (e.g. HOH)
     */
    Atom(const TVector3 v, const double occupancy, const string symbol, const string comp, int serial) {
        // we use our setters so we can validate the input if necessary
        this->set_coordinates(v);
        this->set_occupancy(occupancy);
        this->set_symbol(symbol);
        this->set_comp(comp);
        this->set_serial(serial);
    }

    Atom() {};


    void parse_pdb(const string s) override {
        const char form[] = "ATOM  %5c%4c%1c%3c%1c%4c%1c%8c%8c%8c%6c%6c%2c%2c";
        char serial[5], name[4], altLoc[1], resName[3], chainID[1], resSeq[4], iCode[1], x[8], y[8], z[8], occupancy[6], tempFactor[6], element[2], charge[2];
        sscanf(s.c_str(), form, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge);
        cout << serial << name << altLoc << resName << chainID << resSeq << iCode << x << y << z << occupancy << tempFactor << element << charge << endl;
    }

    /** Calculate the distance to another atom. 
     * @param a the other atom.
     * @return the distance. 
     */
    double distance(const Atom* a) {
        return sqrt(pow(get_x() - a->get_x(), 2) + pow(get_y() - a->get_y(), 2) + pow(get_z() - a->get_z(), 2));
    }

    /** Prints the contents of this object to the terminal. */
    void print() {
        cout << "\nAtom no: " << serial << endl;
        cout << setw(17) << "(x, y, z): (" << setw(6) << get_x() << ", " << setw(6) << get_y() << ", " << setw(6) << get_z() << ")" << endl;
        cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
        cout << setw(16) << "Symbol: " << symbol << endl;
        cout << setw(16) << "Molecule: " << comp << endl;
        return;
    }

    /** Move this atom by a vector.
     * @param v the translation vector.
     */
    void translate(const TVector3 v) {
        coords += v;
    }

    /**
     * @brief Determine if this is a water molecule. 
     * @return true if this is a water molecule, otherwise false. 
     */
    bool is_water() {
        if (comp == "HOH") {
            return true;
        }
        return false;
    }

    // setters
    void set_coordinates(TVector3 v) {this->coords = v;}
    void set_x(double x) {this->coords.SetX(x);}
    void set_y(double y) {this->coords.SetY(y);}
    void set_z(double z) {this->coords.SetZ(z);}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_serial(int serial) {this->serial = serial;}
    void set_comp(string comp) {this->comp = comp;}
    void set_symbol(string symbol) {
        if (!atomic_weight_map.count(symbol)) {
            print_err((format("ERROR: Invalid symbol \"%1%\".") % symbol).str());
            exit(1);
        }
        this->symbol = symbol;
    }

    // getters
    double get_x() const {return coords.X();}
    double get_y() const {return coords.Y();}
    double get_z() const {return coords.Z();}
    TVector3 get_coords() const {return coords;}
    double get_occupancy() const {return occupancy;}
    int get_serial() const {return serial;}
    string get_symbol() const {return symbol;}
    string get_comp() const {return comp;}

    double get_atomic_weight() const {
        if (symbol == "") {
            print_err("ERROR: Attempted to get atomic weight, but the symbol was not set!");
            exit(1);
        }
        return atomic_weight_map.at(symbol);
    };

    /**
     * @brief Comparison function to allow this class to be a map key. 
     * @param rhs Atom to compare against.
     */
    bool operator<(const Atom& rhs) const {
        return serial < rhs.get_serial();
    }

private:
    int serial; // the serial (or sequence number) of this atom. Used when it is read from a file. 

    TVector3 coords; // (x, y, z) coordinates
    double occupancy; // the occupancy (or weight) of this atom
    double atomic_weight; // the atomic weight
    string symbol; // atomic symbol
    string comp = "undefined"; // molecular name, e.g. HOH

    // atomic weights taken from https://www.britannica.com/science/atomic-weight
    const std::map<string, int> atomic_weight_map = {{"H", 1.01}, {"He", 4.00}, {"Li", 6.95}, {"C", 12.01}, {"N", 14.01}, {"O", 16.00}, {"S", 32.06}};
};