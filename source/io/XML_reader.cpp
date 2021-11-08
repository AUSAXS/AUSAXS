#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Reader.h"

class XML_reader : public Reader {
public: 
    /** Constructor for the XML_reader class. 
     * @param filename the name of the input file
     */
    XML_reader(string filename) : Reader(filename) {
        if (filename.find(".xml") == string::npos) {
            print_err("Input file \"" + filename + "\" is not a .xml file!");
            exit(1);
        }
    };

    unique_ptr<File> read() override {
        string line; // placeholder for the current line
        unique_ptr<File> file = std::make_unique<File>();
        while(getline(input, line)) {
            if (line.find("PDBx:atom_site id") == string::npos) {
                file->add("HEADER", line);
                continue; // otherwise we just skip it
            }


            string type = line.substr(0, 6); // read the first 6 characters
            switch(File::get_type(type)) {
                case Record::RecordType::ATOM: {
                    shared_ptr<Atom> atom = std::make_shared<Atom>();
                    atom->parse_xml(line);
                    file->add(atom);
                    break;
                } case Record::RecordType::TERMINATE: {
                    shared_ptr<Terminate> term = std::make_shared<Terminate>();
                    term->parse_xml(line);
                    file->add(term);
                    break;
                } case Record::RecordType::HEADER: {
                    file->add("HEADER", line);
                    break;
                } case Record::RecordType::FOOTER: {
                    file->add("FOOTER", line);
                    break;
                } default: { 
                    print_err("ERROR: Unrecognized type.");
                    exit(1);
                }
            };
        }
        return file;
    } 

    /** Parse the atoms from the input file
     * @return A vector containing all of the parsed atoms.
     */
    vector<Atom*> read2() {
        string line;
        Atom* atom = new Atom();
        vector<Atom*> atoms;
        int serial = 1;
        atom->set_serial(serial);

        while(getline(input, line)) {
            // check if this line contains any relevant information
            if (line.find("PDBx:") == string::npos) {
                continue; // otherwise we just skip it
            }

            // check if the line is the start of an atom_site
            if (line.find("<PDBx:atom_site ") != string::npos) {
                // since this is the start of an atom_site, we make a sanity check on the serial number (it should be the previous serial+1)
                int v_start = line.find("\"")+1;
                int v_end = line.find_last_of("\"");
                string v = line.substr(v_start, v_end-v_start);

                // sanity check
                if (atoi(v.c_str()) != serial) {
                    print_err((format("ERROR: Broken reading sequence in file %1%. Expected id %2%, but found %3%.") % get_filename() % serial % v).str());
                    exit(1);
                }
                continue;
            }

            // check if the line is the end of an atom site
            else if (line.find("</PDBx:atom_site>") != string::npos) {
                // at this point, the atom should contain all of the relevant properties, so we simply store it in an array and prepare another one.
                atoms.push_back(atom);
                atom = new Atom();

                serial++;
                atom->set_serial(serial);
                continue;
            }

            // the remaining options all follow the same format of "<PDBx:(w)>(v)</PDBx:(w)>". we want to find (w) and (v)
            int w_start = line.find("<PDBx:")+6; // the characters in front of w
            int w_end = line.find_first_of(">"); // the character after w
            int v_start = w_end+1; // the start of v
            int v_end = line.find("</PDBx:"); // the characters after v

            string w = line.substr(w_start, w_end-w_start);
            string v = line.substr(v_start, v_end-v_start);

            // check if the word is a coordinate
            if (w == "Cartn_x") {
                atom->set_x(stod(v));
                continue;
            } else if (w == "Cartn_y") {
                atom->set_y(stod(v));
                continue;
            } else if (w == "Cartn_z") {
                atom->set_z(stod(v));
                continue;
            }
            
            // check if the word is the weight
            else if (w == "occupancy") {
                atom->set_occupancy(stod(v));
                continue;
            }

            // check if the word is the symbol
            else if (w == "type_symbol") {
                atom->set_symbol(v);
                continue;
            }

            // check if the word describes the molecule
            else if (w == "label_comp_id") {
                atom->set_comp(v);
                continue;
            }

            // debug
            // cout << "No match found for word: \"" << w << "\" with value \"" << v << "\"" << endl;
        }
        return atoms;
    }
};