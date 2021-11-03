#pragma once

// includes
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>

class Atom {
public:
    /** Constructor for the Atom class. 
     * @param v a vector containing the x, y, z coordinates for the atom. 
     * @param weight the weight of this atom
     * @param symbol the atomic symbol of the base atom
     * @param comp the molecule (e.g. HOH)
     */
    Atom(std::vector<double> v, double weight, std::string symbol, std::string comp) {
        this->set_coordinates(v);
        this->weight = weight;
        this->symbol = symbol;
        this->comp = comp;
    }

    Atom() {};

    /** Calculate the distance to another atom. 
     * @param a the other atom.
     * @return the distance. 
     */
    double distance(Atom* a) {
        return std::sqrt(std::pow(x - a->get_x(), 2) + std::pow(y - a->get_y(), 2) + std::pow(z - a->get_z(), 2));
    }

    /** Returns a standard PDB format representation of this atom.
     * @param serial the entry number of this atom. 
     * @return a PDB string representation of this atom. 
     */
    std::string to_pdb(int serial) {
        return "ATOM  ";
    }

    /** Returns a standard PDBML format representation of this atom. 
     * @return a PDBML string representation of this atom. 
     */
    std::string to_pdbml() {
        return "<PDBx:atom_site id=\"" + std::to_string(serial) + "\"> \
        \n    <PDBx:Cartn_x>" + std::to_string(x) + "</PDBx:Cartn_x> \
        \n    <PDBx:Cartn_y>" + std::to_string(y) + "</PDBx:Cartn_y> \
        \n    <PDBx:Cartn_z>" + std::to_string(z) + "</PDBx:Cartn_z> \
        \n    <PDBx:occupancy>" + std::to_string(weight) + "</PDBx:occupancy> \
        \n    <PDBx:type_symbol>" + symbol + "</PDBx:type_symbol> \
        \n</PDBx:atom_site>";
    }

    /** Prints the contents of this object.
     */
    void print() {
        std::cout << "\nAtom no: " << serial << std::endl;
        std::cout << std::setw(17) << "(x, y, z): (" << std::setw(6) << x << ", " << std::setw(6) << y << ", " << std::setw(6) << z << ")" << std::endl;
        std::cout << std::setw(16) << "Weight: " << std::to_string(weight) << std::endl;
        std::cout << std::setw(16) << "Symbol: " << symbol << std::endl;
        std::cout << std::setw(16) << "Molecule: " << comp << std::endl;
        return;
    }

    // setters
    void set_coordinates(std::vector<double> v) {
        this->x = v[0];
        this->y = v[1];
        this->z = v[2];
    }
    void set_x(double x) {this->x = x;}
    void set_y(double y) {this->y = y;}
    void set_z(double z) {this->z = z;}
    void set_weight(double weight) {this->weight = weight;}
    void set_serial(int serial) {this->serial = serial;}
    void set_comp(std::string comp) {this->comp = comp;}
    void set_symbol(std::string symbol) {this->symbol = symbol;}

    // getters
    int get_x() {return x;}
    int get_y() {return y;}
    int get_z() {return z;}
    int get_weight() {return weight;}
    int get_serial() {return serial;}
    std::string get_symbol() {return symbol;}
    std::string get_comp() {return comp;}

private:
    int serial; // the serial (or sequence number) of this atom. Used when it is read from a file. 
    double x, y, z; // coordinates
    double weight; // weight
    std::string symbol; // atomic symbol
    std::string comp = "undefined"; // special name e.g. Ca

    const std::map<std::string, int> atomic_number_map = {{"H", 1}, {"He", 2}, {"C", 6}};
};