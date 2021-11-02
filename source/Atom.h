#pragma once

// includes
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>

class Atom {
public:
    int serial; // the serial (or sequence number) of this atom. Used when it is read from a file. 
    double x, y, z; // coordinates
    double w; // weight
    std::string symbol; // atomic symbol
    std::string name; // special name e.g. Ca

    /** Constructor for the Atom class. 
     * @param v a vector containing the x, y, z coordinates for the atom. 
     * @param w the weight of this atom
     * @param symbol the atomic symbol of the base atom
     * @param name the special? name
     */
    Atom(std::vector<double> v, double w, std::string symbol, std::string name) {
        this->set_coordinates(v);
        this->w = w;
        this->symbol = symbol;
        this->name = name;
    }

    Atom() {};

    /** Return a copy of the coordinates of this atom. 
     * @return a vector containing the x, y, z coordinates. 
     */
    std::vector<double> get_coordinates() {
        std::vector<double> v = {x, y, z};
        return v;
    }

    /** Set the coordinates for this atom. 
     * @param v The x, y, z coordinates for the new location. 
     */
    void set_coordinates(std::vector<double> v) {
        this->x = v[0];
        this->y = v[1];
        this->z = v[2];
    }

    /** Set the x-coordinate. 
     * @param x the new x-coordinate.
     */
    void set_x(double x) {
        this->x = x;
    }

    /** Set the y-coordinate. 
     * @param x the new y-coordinate.
     */
    void set_y(double y) {
        this->y = y;
    }

    /** Set the z-coordinate. 
     * @param x the new z-coordinate.
     */
    void set_z(double z) {
        this->z = z;
    }

    /** Set the weight. 
     * @param w the new x-coordinate.
     */
    void set_w(double w) {
        this->w = w;
    }

    /** Set the symbol. 
     * @param symbol the new symbol.
     */
    void set_symbol(std::string symbol) {
        this->symbol = symbol;
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
        \n    <PDBx:occupancy>" + std::to_string(w) + "</PDBx:occupancy> \
        \n    <PDBx:type_symbol>" + symbol + "</PDBx:type_symbol> \
        \n</PDBx:atom_site>";
    }

    /** Prints the contents of this object.
     */
    void print() {
        std::cout << "\nAtom no: " << serial << std::endl;
        std::cout << std::setw(17) << "(x, y, z): (" << std::setw(6) << x << ", " << std::setw(6) << y << ", " << std::setw(6) << z << ")" << std::endl;
        std::cout << std::setw(16) << "Weight: " << std::to_string(w) << std::endl;
        std::cout << std::setw(16) << "Symbol: " << symbol << std::endl;
        return;
    }

private:
    const std::map<std::string, int> atomic_number_map = {{"H", 1}, {"He", 2}, {"C", 6}};
};