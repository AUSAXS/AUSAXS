#include <io/ProteinFile.h>
#include <io/ExistingFile.h>
#include <io/PDBReader.h>
#include <io/PDBWriter.h>
#include <data/Record.h>
#include <data/Protein.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <utility/Exceptions.h>

ProteinFile::ProteinFile() = default;

ProteinFile::ProteinFile(const ProteinFile& ProteinFile) : header(ProteinFile.header), footer(ProteinFile.footer), terminate(ProteinFile.terminate), 
    protein_atoms(ProteinFile.protein_atoms), hydration_atoms(ProteinFile.hydration_atoms) {}

ProteinFile::ProteinFile(ProteinFile&& ProteinFile) noexcept : header(ProteinFile.header), footer(ProteinFile.footer), terminate(ProteinFile.terminate), 
    protein_atoms(std::move(ProteinFile.protein_atoms)), hydration_atoms(std::move(ProteinFile.hydration_atoms)) {}

ProteinFile::ProteinFile(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms) : protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

ProteinFile::ProteinFile(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms, const Header& header, const Footer& footer, const Terminate& terminate) 
    : header(header), footer(footer), terminate(terminate), protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

ProteinFile::ProteinFile(const io::ExistingFile& ProteinFilename) {
    reader = construct_reader(ProteinFilename);
    read(ProteinFilename);
}

ProteinFile::~ProteinFile() = default;

void ProteinFile::update(std::vector<Atom>& patoms, std::vector<Water>& hatoms) {
    protein_atoms = patoms;
    hydration_atoms = hatoms;
}

void ProteinFile::read(const io::ExistingFile& path) {reader->read(path);}

void ProteinFile::write(const io::File& path) {
    refresh();
    writer = construct_writer(path);
    writer->write(path);
}

const std::vector<Atom>& ProteinFile::get_protein_atoms() const {return protein_atoms;}

const std::vector<Water> ProteinFile::get_hydration_atoms() const {return hydration_atoms;}

void ProteinFile::add(const Atom& r) {
    protein_atoms.push_back(r);
}

void ProteinFile::add(Atom&& r) {
    protein_atoms.push_back(std::move(r));
}

void ProteinFile::add(const Water& r) {
    hydration_atoms.push_back(r);
}

void ProteinFile::add(Water&& r) {
    hydration_atoms.push_back(std::move(r));
}

void ProteinFile::add(const Terminate& ter) {
    terminate = ter;
    // terminates.push_back(r);
}

void ProteinFile::add(const RecordType& type, const std::string& s) {
    if (type == RecordType::HEADER) {
        header.add(s);
    } else if (type == RecordType::FOOTER) {
        footer.add(s);
    } else {
        throw except::invalid_argument("ProteinFile::add: Type is not \"HEADER\" or \"FOOTER\"!");
    }
}

ProteinFile ProteinFile::copy() const {
    return ProteinFile(*this);
}

ProteinFile& ProteinFile::operator=(const ProteinFile& rhs) {
    protein_atoms = rhs.protein_atoms;
    hydration_atoms = rhs.hydration_atoms;
    header = rhs.header;
    footer = rhs.footer;
    terminate = rhs.terminate;
    return *this;
}

ProteinFile& ProteinFile::operator=(ProteinFile&& rhs) = default;

void ProteinFile::refresh() {
    if (protein_atoms.empty()) {return;}

    bool terminate_inserted = false;
    char chainID = '0'; int resSeq = 0; int serial = protein_atoms[0].serial;

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
        if (!terminate_inserted && a.get_type() == RecordType::WATER) {
            insert_ter();
            resSeq++; // TER records always denotes the end of a sequence
            serial++;
        }
        a.set_serial(serial++ % 100000); // fix possible errors in the serial
    }

    if (!terminate_inserted) {
        insert_ter();
        resSeq++; // TER records always denotes the end of a sequence
        serial++;
    }

    chainID = protein_atoms[protein_atoms.size()-1].chainID+1;
    resSeq = protein_atoms[protein_atoms.size()-1].resSeq + 1;
    for (auto& a : hydration_atoms) {
        a.set_serial(serial++ % 100000);
        a.set_resSeq(resSeq++ % 10000);
        a.set_chainID(chainID + int(resSeq/10000));
    }
}

bool ProteinFile::operator==(const ProteinFile& rhs) const = default;

#define FAILURE_MSG false
#if FAILURE_MSG
    #include <iostream>
#endif
bool ProteinFile::equals_content(const ProteinFile& rhs) const {
    if (protein_atoms.size() != rhs.protein_atoms.size()) {
        #if FAILURE_MSG
            std::cout << "protein_atoms.size() != rhs.protein_atoms.size()" << std::endl;
        #endif
        return false;
    }

    if (hydration_atoms.size() != rhs.hydration_atoms.size()) {
        #if FAILURE_MSG
            std::cout << "hydration_atoms.size() != rhs.hydration_atoms.size()" << std::endl;
        #endif
        return false;
    }

    for (unsigned int i = 0; i < protein_atoms.size(); i++) {
        if (!protein_atoms[i].equals_content(rhs.protein_atoms[i])) {
            #if FAILURE_MSG
                std::cout << "!protein_atoms[" << i << "].equals_content(rhs.protein_atoms[" << i << "])" << std::endl;
            #endif
            return false;
        }
    }

    for (unsigned int i = 0; i < hydration_atoms.size(); i++) {
        if (!hydration_atoms[i].equals_content(rhs.hydration_atoms[i])) {
            #if FAILURE_MSG
                std::cout << "!hydration_atoms[" << i << "].equals_content(rhs.hydration_atoms[" << i << "])" << std::endl;
            #endif
            return false;
        }
    }

    return true;    
}

std::unique_ptr<Reader> ProteinFile::construct_reader(const io::ExistingFile& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml ProteinFile
        throw except::invalid_argument("ProteinFile::construct_reader: .xml input ProteinFiles are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb ProteinFile
        return std::make_unique<PDBReader>(this);
    } else if (path.extension() == ".ent" || path.extension() == ".ENT") { // .ent ProteinFile
        return std::make_unique<PDBReader>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("ProteinFile::construct_reader: Unsupported ProteinFile extension of input ProteinFile \"" + path + "\".");
    }
}

std::unique_ptr<Writer> ProteinFile::construct_writer(const io::File& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml ProteinFile
        throw except::invalid_argument("ProteinFile::construct_writer: .xml input ProteinFiles are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb ProteinFile
        return std::make_unique<PDBWriter>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("ProteinFile::construct_writer: Unsupported ProteinFile extension of input ProteinFile \"" + path + "\".");
    }
}