#pragma once

// includes
#include <string>
#include <vector>
#include <utility>
#include <fstream>

// my own includes
#include "data/Record.h"
#include "data/Terminate.h"
#include "data/Header.h"
#include "data/Footer.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "data/File.h"
#include "settings.h"

using std::vector, std::string, std::cout, std::endl, std::unique_ptr, std::shared_ptr; 

class PDB_file : public File {
public: 
    /** 
     * @brief Constructor for the PDB_file class. 
     */
    PDB_file(string filename) : File(filename) {
        read();
    }
    
    /**
     * @brief write this File to disk. 
     * @param path the output path.
     */
    void write(string path) const override {
        std::ofstream output(path);
        if (!output.is_open()) {
            print_err("Error in PDB_file::write: Could not open file \"" + path + "\"");
            exit(1);
        }
        output << as_pdb() << std::flush;
        output.close();
        cout << "Output written to file " + path + "." << endl;
    }

    /**
     * @brief Update the contents of this File to reflect the input Protein.
     * @param protein the protein to use.
     */
    void update(vector<shared_ptr<Atom>> protein_atoms, vector<shared_ptr<Hetatom>> hydration_atoms) override {
        vector<shared_ptr<Record>> c(protein_atoms.size() + hydration_atoms.size() + 1, nullptr);
        int i = 0; // counter

        // insert all ATOMs
        bool terminate_inserted = false;
        string chainID = "0"; int resSeq = 0;
        auto insert_ter = [&] () {
            // last atom before the terminate
            // we need this to determine what chainID and resSeq to use for the terminate and hetatms
            shared_ptr<Atom> a = std::static_pointer_cast<Atom>(c[i-1]);
            chainID = a->get_chainID();
            resSeq = a->get_resSeq();

            // insert a TER record between the atoms and waters
            if (i != 0) {
                c[i] = std::make_shared<Terminate>(i+1, a->get_resName(), chainID, resSeq, " ");
            }
            terminate_inserted = true;
        };

        for (auto const& a : protein_atoms) {
            if (!terminate_inserted && a->get_type() == Record::RecordType::HETATM) {
                insert_ter();
                i++;
            }
            a->set_serial(i+1); // fix possible errors in the serial
            c[i] = a;
            i++;
        }

        if (!terminate_inserted) {
            insert_ter();
            i++;
        }

        // insert all waters
        for (auto const& a : hydration_atoms) {
            a->set_serial(i+1); // fix possible errors in the serial
            a->set_resSeq(resSeq+1);
            a->set_chainID(chainID);
            c[i] = a;
            resSeq++;
            i++;
        }

        contents = c;
    };

private:
    /**
     * @brief Read the file backing this File object. 
     */
    void read() override {
        // check if file was succesfully opened
        std::ifstream input(filename);
        if (!input.is_open()) {
            print_err("Error in PDB_file::read: Could not open file \"" + filename + "\"");
            exit(1);
        }

        string line; // placeholder for the current line
        while(getline(input, line)) {
            string type = line.substr(0, 6); // read the first 6 characters
            switch(Record::get_type(type)) {
                case Record::RecordType::HETATM: {
                    shared_ptr<Hetatom> atom = std::make_shared<Hetatom>();
                    atom->parse_pdb(line);
                    add(atom);
                    break;
                } case Record::RecordType::ATOM: {
                    shared_ptr<Atom> atom = std::make_shared<Atom>();
                    atom->parse_pdb(line);
                    add(atom);
                    break;
                } case Record::RecordType::TERMINATE: {
                    shared_ptr<Terminate> term = std::make_shared<Terminate>();
                    term->parse_pdb(line);
                    add(term);
                    break;
                } case Record::RecordType::HEADER: {
                    add("HEADER", line);
                    break;
                } case Record::RecordType::FOOTER: {
                    add("FOOTER", line);
                    break;
                } default: { 
                    print_err("Error in PDB_file::read: Unrecognized type.");
                    exit(1);
                }
            };
        }
        input.close();
    }
};