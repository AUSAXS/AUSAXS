#pragma once

// includes
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <boost/format.hpp>

// ROOT
#include <TVector3.h>

// my own stuff
#include "Tools.cpp"

using namespace ROOT;
using std::vector, std::string, std::cout, std::endl, std::setw;
using boost::format;

class Atom {
public:
    /** Constructor for the Atom class. 
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom
     * @param symbol the atomic symbol of the base atom
     * @param comp the molecule (e.g. HOH)
     */
    Atom(TVector3 v, double occupancy, string symbol, string comp) {
        this->set_coordinates(v);
        this->occupancy = occupancy;
        this->symbol = symbol;
        this->comp = comp;
    }

    Atom() {};

    /** Calculate the distance to another atom. 
     * @param a the other atom.
     * @return the distance. 
     */
    double distance(Atom* a) {
        return sqrt(pow(get_x() - a->get_x(), 2) + pow(get_y() - a->get_y(), 2) + pow(get_z() - a->get_z(), 2));
    }

    /** Returns a standard PDB format representation of this atom.
     * @param serial the entry number of this atom. 
     * @return a PDB string representation of this atom. 
     */
    string to_pdb(int serial) {
        return "ATOM  ";
    }

    /** Returns a standard PDBML format representation of this atom. 
     * @return a PDBML string representation of this atom. 
     */
    string to_pdbml() {
        return (format("<PDBx:atom_site id=\"%1%\"> \
        \n    <PDBx:Cartn_x>%2%</PDBx:Cartn_x> \
        \n    <PDBx:Cartn_y>%3%</PDBx:Cartn_y> \
        \n    <PDBx:Cartn_z>%4%</PDBx:Cartn_z> \
        \n    <PDBx:occupancy>%5%</PDBx:occupancy> \
        \n    <PDBx:type_symbol>%6%</PDBx:type_symbol> \
        \n</PDBx:atom_site>") % serial % get_x() % get_y() % get_z() % occupancy % symbol).str();
    }

    /** Prints the contents of this object. */
    void print() {
        cout << "\nAtom no: " << serial << endl;
        cout << setw(17) << "(x, y, z): (" << setw(6) << get_x() << ", " << setw(6) << get_y() << ", " << setw(6) << get_z() << ")" << endl;
        cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
        cout << setw(16) << "Symbol: " << symbol << endl;
        cout << setw(16) << "Molecule: " << comp << endl;
        return;
    }

    /** Move this atom by a vector
     * @param v the translation vector.
     */
    void translate(TVector3 v) {
        coords += v;
    }

    // setters
    void set_coordinates(TVector3 v) {this->coords = v;}
    void set_x(double x) {this->coords.SetX(x);}
    void set_y(double y) {this->coords.SetY(y);}
    void set_z(double z) {this->coords.SetZ(z);}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_serial(int serial) {this->serial = serial;}
    void set_comp(string comp) {this->comp = comp;}
    void set_symbol(string symbol) {this->symbol = symbol;}

    // getters
    int get_x() {return coords.X();}
    int get_y() {return coords.Y();}
    int get_z() {return coords.Z();}
    TVector3 get_coords() {return coords;}
    int get_occupancy() {return occupancy;}
    int get_serial() {return serial;}
    string get_symbol() {return symbol;}
    string get_comp() {return comp;}

    double get_atomic_weight() {
        if (symbol == "") {
            print_err("ERROR: Attempted to get atomic weight, but the symbol was not set!");
            exit(1);
        }
        if (!atomic_weight_map.count(symbol)) {
            print_err((format("ERROR: Weight is undefined for the atom \"%1%\"") % symbol).str());
            exit(1);
        }
        return atomic_weight_map.at(symbol);
    };

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