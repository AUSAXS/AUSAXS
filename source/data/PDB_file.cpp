#pragma once

// includes
#include <string>
#include <vector>
#include <utility>

// my own includes
#include "data/Record.h"
#include "data/Terminate.cpp"
#include "data/Header.cpp"
#include "data/Footer.cpp"
#include "data/Atom.cpp"
#include "data/File.h"

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
    virtual void update(vector<shared_ptr<Atom>> protein_atoms, vector<shared_ptr<Atom>> hydration_atoms) {
        vector<shared_ptr<Record>> c(protein_atoms.size() + hydration_atoms.size() + 1, nullptr);
        int i = 0; // counter

        // insert all ATOMs
        for (auto const& a : protein_atoms) {
            a->set_serial(i+1); // fix possible errors in the serial
            c[i] = a;
            i++;
        }
        // last atom before the terminate
        // we need this to determine what chainID and resSeq to use for the terminate and hetatms
        shared_ptr<Atom> a = std::static_pointer_cast<Atom>(c[i-1]);
        string chainID = a->get_chainID();
        int resSeq = a->get_resSeq();

        // insert a TER record between the ATOMs and HETATMs
        if (i != 0) {
            c[i] = find_ter_separator(); // see if we can reuse the original
            if (c[i] == nullptr) { // it could not be found, so we create a new one
                c[i] = std::make_shared<Terminate>(i+1, a->get_resName(), chainID, resSeq, " ");
            } else { // it does exist, so we update its serial and reuse it
                std::static_pointer_cast<Terminate>(c[i])->set_serial(i+1);
            }
            i++;
        }

        // insert all HETATMs
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
     * @brief Find the Terminate Record separating ATOMs and HETATMs
     * @return A pointer to the Terminate Record if it exists, nullptr otherwise.
     */
    shared_ptr<Record> find_ter_separator() {
        for (auto const& r : contents) {
            if (r->get_type() == Record::TERMINATE) {
                return r;
            }
        }
        return nullptr;
    }

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
                case Record::RecordType::ATOM: {
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