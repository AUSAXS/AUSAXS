#include <data/detail/AtomCollection.h>
#include <data/record/Record.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/record/Terminate.h>
#include <data/record/Header.h>
#include <data/record/Footer.h>
#include <data/Molecule.h>
#include <io/ExistingFile.h>
#include <io/PDBReader.h>
#include <io/PDBWriter.h>
#include <utility/Exceptions.h>

using namespace data::detail;

AtomCollection::AtomCollection() = default;

AtomCollection::AtomCollection(const AtomCollection& AtomCollection) : header(AtomCollection.header), footer(AtomCollection.footer), terminate(AtomCollection.terminate), 
    protein_atoms(AtomCollection.protein_atoms), hydration_atoms(AtomCollection.hydration_atoms) {}

AtomCollection::AtomCollection(AtomCollection&& AtomCollection) noexcept : header(AtomCollection.header), footer(AtomCollection.footer), terminate(AtomCollection.terminate), 
    protein_atoms(std::move(AtomCollection.protein_atoms)), hydration_atoms(std::move(AtomCollection.hydration_atoms)) {}

AtomCollection::AtomCollection(const std::vector<record::Atom>& protein_atoms, const std::vector<record::Water>& hydration_atoms) : protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

AtomCollection::AtomCollection(const std::vector<record::Atom>& protein_atoms, const std::vector<record::Water>& hydration_atoms, const record::Header& header, const record::Footer& footer, const record::Terminate& terminate) 
    : header(header), footer(footer), terminate(terminate), protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

AtomCollection::AtomCollection(const io::File& file) {
    reader = construct_reader(file);
    read(file);
}

AtomCollection::~AtomCollection() = default;

void AtomCollection::update(std::vector<record::Atom>& patoms, std::vector<record::Water>& hatoms) {
    protein_atoms = patoms;
    hydration_atoms = hatoms;
}

void AtomCollection::read(const io::File& path) {reader->read(path);}

void AtomCollection::write(const io::File& path) {
    refresh();
    writer = construct_writer(path);
    writer->write(path);
}

const std::vector<data::record::Atom>& AtomCollection::get_protein_atoms() const {return protein_atoms;}

const std::vector<data::record::Water> AtomCollection::get_hydration_atoms() const {return hydration_atoms;}

void AtomCollection::add(const record::Atom& r) {
    protein_atoms.push_back(r);
}

void AtomCollection::add(record::Atom&& r) {
    protein_atoms.push_back(std::move(r));
}

void AtomCollection::add(const record::Water& r) {
    hydration_atoms.push_back(r);
}

void AtomCollection::add(record::Water&& r) {
    hydration_atoms.push_back(std::move(r));
}

void AtomCollection::add(const record::Terminate& ter) {
    terminate = ter;
    // terminates.push_back(r);
}

void AtomCollection::add(const record::RecordType& type, const std::string& s) {
    if (type == record::RecordType::HEADER) {
        header.add(s);
    } else if (type == record::RecordType::FOOTER) {
        footer.add(s);
    } else {
        throw except::invalid_argument("AtomCollection::add: Type is not \"HEADER\" or \"FOOTER\"!");
    }
}

AtomCollection AtomCollection::copy() const {
    return AtomCollection(*this);
}

AtomCollection& AtomCollection::operator=(const AtomCollection& rhs) {
    protein_atoms = rhs.protein_atoms;
    hydration_atoms = rhs.hydration_atoms;
    header = rhs.header;
    footer = rhs.footer;
    terminate = rhs.terminate;
    return *this;
}

AtomCollection& AtomCollection::operator=(AtomCollection&& rhs) = default;

void AtomCollection::refresh() {
    if (protein_atoms.empty()) {return;}

    bool terminate_inserted = false;
    char chainID = '0'; int resSeq = 0; int serial = protein_atoms[0].serial;

    auto insert_ter = [&] () {
        // last atom before the terminate
        // we need this to determine what chainID and resSeq to use for the terminate and hetatms
        const record::Atom& a = protein_atoms.at(serial-1-protein_atoms[0].serial);
        chainID = a.chainID;
        resSeq = a.resSeq;
        if (serial != 0) {terminate = record::Terminate(serial, a.resName, a.chainID, a.resSeq, " ");}
        terminate_inserted = true;
    };

    for (auto& a : protein_atoms) {
        if (!terminate_inserted && a.get_type() == record::RecordType::WATER) {
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
        a.set_residue_sequence_number(resSeq++ % 10000);
        a.set_chainID(chainID + int(resSeq/10000));
    }
}

bool AtomCollection::operator==(const AtomCollection& rhs) const = default;

#define FAILURE_MSG false
#if FAILURE_MSG
    #include <iostream>
#endif
bool AtomCollection::equals_content(const AtomCollection& rhs) const {
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

std::unique_ptr<io::detail::Reader> AtomCollection::construct_reader(const io::File& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml AtomCollection
        throw except::invalid_argument("AtomCollection::construct_reader: .xml input AtomCollections are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb AtomCollection
        return std::make_unique<io::detail::PDBReader>(this);
    } else if (path.extension() == ".ent" || path.extension() == ".ENT") { // .ent AtomCollection
        return std::make_unique<io::detail::PDBReader>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("AtomCollection::construct_reader: Unsupported AtomCollection extension of input AtomCollection \"" + path + "\".");
    }
}

std::unique_ptr<io::detail::Writer> AtomCollection::construct_writer(const io::File& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml AtomCollection
        throw except::invalid_argument("AtomCollection::construct_writer: .xml input AtomCollections are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb AtomCollection
        return std::make_unique<io::detail::PDBWriter>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("AtomCollection::construct_writer: Unsupported AtomCollection extension of input AtomCollection \"" + path + "\".");
    }
}