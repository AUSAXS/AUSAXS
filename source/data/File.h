#pragma once

// includes
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

// my own includes
#include "data/Record.h"
#include "data/Terminate.h"
#include "data/Header.h"
#include "data/Footer.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

using std::vector, std::string, std::unique_ptr, std::shared_ptr; 

class File {
public: 
    File(string filename) : filename(filename) {}
    virtual ~File() {}

    Header header;
    Footer footer;
    Terminate terminate;
    vector<Atom> protein_atoms;
    vector<Hetatom> hydration_atoms;

    /**
     * @brief Update the contents of this File to reflect the input Protein.
     * @param protein the protein to use.
     */
    void update(vector<Atom>& patoms, vector<Hetatom>& hatoms) {
        protein_atoms = patoms;
        hydration_atoms = hatoms;
    }

    /**
     * @brief Read the file backing this File object. 
     */
    virtual void read() = 0;

    /**
     * @brief write this File to disk. 
     * @param path the output path.
     */
    virtual void write(const string path) = 0;

    /**
     * @brief Get the protein atoms contained in this File. 
     */
    const vector<Atom>& get_protein_atoms() const {return protein_atoms;}

    /**
     * @brief Get the hydration atoms contained in this File. 
     */
    const vector<Hetatom> get_hydration_atoms() const {return hydration_atoms;}

    /**
     * @brief Create a string representation of this File.
     * @return The string representation. 
     */
    string as_pdb() const {
        string s;
        s += header.get();

        size_t i_ter = terminate.serial;
        bool printed_ter = false;
        for (size_t i = 0; i < protein_atoms.size(); i++) {
            if (i == i_ter) { // check if this is where the terminate is supposed to go
                s += terminate.as_pdb(); // write it if so
                printed_ter = true;
            }
            s += protein_atoms[i].as_pdb();
        }
        if (!printed_ter) {s += terminate.as_pdb();}
        for (size_t i = 0; i < hydration_atoms.size(); i++) {s += hydration_atoms[i].as_pdb();}

        s += footer.get();
        return s;
    };

    /** 
     * @brief Add an Atom record to this File. 
     * @param r Atom to be added. 
     */
    virtual void add(const Atom r) {
        protein_atoms.push_back(r);
    }

    /** 
     * @brief Add a Hetatom record to this File. 
     * @param r Hetatom to be added. 
     */
    virtual void add(const Hetatom r) {
        hydration_atoms.push_back(r);
    }

    /**
     * @brief Add a Terminate record to this File. 
     * @param r Terminate to be added. 
     */
    void add(const Terminate) {
        // terminates.push_back(r);
    }

    /**
     * @brief Add a header or footer record to this File. 
     * @param type HEADER or FOOTER.
     * @param s string to be added. 
     */
    void add(const string type, const string s) {
        if (type == "HEADER") {
            header.add(s);
        } else if (type == "FOOTER") {
            footer.add(s);
        } else {
            print_err("Error in File::add: string " + type + " is not \"HEADER\" or \"FOOTER\"!");
            exit(1);
        }
    }

protected:
    /**
     * @brief This method guarantees that the File is in a valid state for printing.
     */
    void refresh() {
        bool terminate_inserted = false;
        string chainID = "0"; int resSeq = 0; int serial = protein_atoms[0].get_serial();

        auto insert_ter = [&] () {
            // last atom before the terminate
            // we need this to determine what chainID and resSeq to use for the terminate and hetatms
            const Atom& a = protein_atoms.at(serial-1-protein_atoms[0].get_serial());
            chainID = a.get_chainID();
            resSeq = a.get_resSeq();
            if (serial != 0) {terminate = Terminate(serial, a.get_resName(), a.get_chainID(), a.get_resSeq(), " ");}
            terminate_inserted = true;
        };

        for (auto& a : protein_atoms) {
            if (!terminate_inserted && a.get_type() == Record::RecordType::HETATM) {
                insert_ter();
                resSeq++; // TER records always denotes the end of a sequence
                serial++;
            }
            a.set_serial(serial); // fix possible errors in the serial
            serial++;
        }

        if (!terminate_inserted) {
            insert_ter();
            resSeq++; // TER records always denotes the end of a sequence
            serial++;
        }

        chainID = protein_atoms[protein_atoms.size()-1].get_chainID();
        resSeq = protein_atoms[protein_atoms.size()-1].get_resSeq() + 1;
        for (auto& a : hydration_atoms) {
            a.set_serial(serial);
            a.set_resSeq(resSeq);
            a.set_chainID(chainID);
            resSeq++;
            serial++;
        }
    }

    string filename;
};