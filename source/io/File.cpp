#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include "io/File.h"

File::File(const File& file) : header(file.header), footer(file.footer), terminate(file.terminate), 
    protein_atoms(file.protein_atoms), hydration_atoms(file.hydration_atoms) {}

File::File(const File&& file) noexcept : header(file.header), footer(file.footer), terminate(file.terminate), 
    protein_atoms(std::move(file.protein_atoms)), hydration_atoms(std::move(file.hydration_atoms)) {}

File::File(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms) : protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

File::File(string filename) {
    reader = construct_reader(filename);
    read(filename);
}

File::~File() = default;

void File::update(vector<Atom>& patoms, vector<Hetatom>& hatoms) {
    protein_atoms = patoms;
    hydration_atoms = hatoms;
}

void File::read(const string& path) {reader->read(path);}

void File::write(const string path) {
    writer = construct_writer(path);
    writer->write(path);
}

const vector<Atom>& File::get_protein_atoms() const {return protein_atoms;}

const vector<Hetatom> File::get_hydration_atoms() const {return hydration_atoms;}

void File::add(const Atom r) {
    protein_atoms.push_back(r);
}

void File::add(const Hetatom r) {
    hydration_atoms.push_back(r);
}

void File::add(const Terminate) {
    // terminates.push_back(r);
}

void File::add(const string type, const string s) {
    if (type == "HEADER") {
        header.add(s);
    } else if (type == "FOOTER") {
        footer.add(s);
    } else {
        print_err("Error in File::add: string " + type + " is not \"HEADER\" or \"FOOTER\"!");
        exit(1);
    }
}

File File::copy() const {
    return File(*this);
}

void File::refresh() {
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

std::unique_ptr<Reader> File::construct_reader(const string& path) {
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

std::unique_ptr<Writer> File::construct_writer(const string& path) {
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