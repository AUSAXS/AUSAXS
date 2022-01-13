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
#include "io/Reader.h"
#include "io/Writer.h"
#include "io/PDBWriter.h"
#include "io/PDBReader.h"

using std::vector, std::string, std::unique_ptr, std::shared_ptr; 

/**
 * @brief \class File
 * An abstract representation of an input file. Handles both reading and writing of files, and also stores the relevant file data. 
 */
class File {
public: 
    /**
     * @brief Constructor. 
     *        Construct a new File based on two vectors of atoms. 
     * @param protein_atoms A vector containing the constituent atoms of the molecule. 
     * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
     */
    File(const vector<Atom> protein_atoms, const vector<Hetatom> hydration_atoms) : protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

    /**
     * @brief Constructor.
     *        Construct a new File based on a input molecular data file. 
     * @param filename Path to the input file. 
     */
    File(string filename) {
        reader = construct_reader(filename);
        read(filename);
    }

    /**
     * @brief Destructor.
     */
    virtual ~File() {}

    /**
     * @brief Update the contents of this file.
     * @param patoms The new constituent atoms of the molecule. 
     * @param hatoms The new hydration layer. 
     */
    void update(vector<Atom>& patoms, vector<Hetatom>& hatoms) {
        protein_atoms = patoms;
        hydration_atoms = hatoms;
    }

    /**
     * @brief Fill this object with data from a given input data file. 
     * @param path Path to the input data file. 
     */
    virtual void read(const string& path) {reader->read(path);}

    /**
     * @brief Write this file to disk. 
     * @param path Path to where this object will be written. 
     */
    virtual void write(const string path) {
        writer = construct_writer(path);
        writer->write(path);
    }

    /**
     * @brief Get the protein atoms contained in this File. 
     */
    const vector<Atom>& get_protein_atoms() const {return protein_atoms;}

    /**
     * @brief Get the hydration atoms contained in this File. 
     */
    const vector<Hetatom> get_hydration_atoms() const {return hydration_atoms;}

    /** 
     * @brief Add an Atom record to this file. 
     * @param r The Atom record to be added.
     */
    virtual void add(const Atom r) {
        protein_atoms.push_back(r);
    }

    /** 
     * @brief Add a Hetatom record to this file. 
     * @param r The Hetatom record to be added. 
     */
    virtual void add(const Hetatom r) {
        hydration_atoms.push_back(r);
    }

    /**
     * @brief Add a Terminate record to this file. 
     * @param r The Terminate record to be added. 
     */
    void add(const Terminate) {
        // terminates.push_back(r);
    }

    /**
     * @brief Add a header or footer record to this file. 
     * @param type HEADER or FOOTER.
     * @param s The text string to be added. 
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

    /**
     * @brief Internally updates the consistency of the currently stored data. This ensures this object is in a valid
     *        state for printing. 
     */
    void refresh() {
        bool terminate_inserted = false;
        string chainID = "0"; int resSeq = 0; int serial = protein_atoms[0].serial;

        auto insert_ter = [&] () {
            // last atom before the terminate
            // we need this to determine what chainID and resSeq to use for the terminate and hetatms
            const Atom& a = protein_atoms.at(serial-1-protein_atoms[0].serial);
            chainID = a.chainID;
            resSeq = a.resSeq;
            if (serial != 0) {terminate = Terminate(serial, a.resName, a.chainID, a.resSeq, " ");}
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

        chainID = protein_atoms[protein_atoms.size()-1].chainID;
        resSeq = protein_atoms[protein_atoms.size()-1].resSeq + 1;
        for (auto& a : hydration_atoms) {
            a.set_serial(serial);
            a.set_resSeq(resSeq);
            a.set_chainID(chainID);
            resSeq++;
            serial++;
        }
    }

    Header header;
    Footer footer;
    Terminate terminate;
    vector<Atom> protein_atoms;
    vector<Hetatom> hydration_atoms;

private:
    std::unique_ptr<Reader> reader;
    std::unique_ptr<Writer> writer;

    /**
     * @brief Construct a Reader appropriate for the file format deduced from the input data file. 
     * @param path Path to the input data file. 
     */
    std::unique_ptr<Reader> construct_reader(const string& path) {
        if (path.find(".xml") != string::npos || path.find(".XML") != string::npos) { // .xml file
            print_err("Error in Protein::Protein: .xml input files are not supported.");
            exit(1);
        } else if (path.find(".pdb") != string::npos || path.find(".PDB") != string::npos) { // .pdb file
            return std::make_unique<PDBReader>(this);
        } else { // anything else - we cannot handle this
            print_err((format("Error in Protein::Protein: Invalid file extension of input file %1%.") % path).str());
            exit(1);
        }
    }

    /**
     * @brief Construct a Writer appropriate for the file format deduced from the output save path.
     * @param path Path to where this object will be written. 
     */
    std::unique_ptr<Writer> construct_writer(const string& path) {
        if (path.find(".xml") != string::npos || path.find(".XML") != string::npos) { // .xml file
            print_err("Error in Protein::Protein: .xml input files are not supported.");
            exit(1);
        } else if (path.find(".pdb") != string::npos || path.find(".PDB") != string::npos) { // .pdb file
            return std::make_unique<PDBWriter>(this);
        } else { // anything else - we cannot handle this
            print_err((format("Error in Protein::Protein: Invalid file extension of input file %1%.") % path).str());
            exit(1);
        }
    }
};