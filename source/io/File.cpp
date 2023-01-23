#include <io/File.h>

File::File(const File& file) : header(file.header), footer(file.footer), terminate(file.terminate), 
    protein_atoms(file.protein_atoms), hydration_atoms(file.hydration_atoms) {}

File::File(File&& file) noexcept : header(file.header), footer(file.footer), terminate(file.terminate), 
    protein_atoms(std::move(file.protein_atoms)), hydration_atoms(std::move(file.hydration_atoms)) {}

File::File(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms) : protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

File::File(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms, const Header& header, const Footer& footer, const Terminate& terminate) 
    : header(header), footer(footer), terminate(terminate), protein_atoms(protein_atoms), hydration_atoms(hydration_atoms) {}

File::File(std::string filename) {
    reader = construct_reader(filename);
    read(filename);
}

File::~File() = default;

void File::update(std::vector<Atom>& patoms, std::vector<Water>& hatoms) {
    protein_atoms = patoms;
    hydration_atoms = hatoms;
}

void File::read(std::string path) {reader->read(path);}

void File::write(std::string path) {
    refresh();
    writer = construct_writer(path);
    writer->write(path);
}

const std::vector<Atom>& File::get_protein_atoms() const {return protein_atoms;}

const std::vector<Water> File::get_hydration_atoms() const {return hydration_atoms;}

void File::add(const Atom& r) {
    protein_atoms.push_back(r);
}

void File::add(const Water& r) {
    hydration_atoms.push_back(r);
}

void File::add(const Terminate& ter) {
    terminate = ter;
    // terminates.push_back(r);
}

void File::add(Record::RecordType type, std::string s) {
    if (type == Record::RecordType::HEADER) {
        header.add(s);
    } else if (type == Record::RecordType::FOOTER) {
        footer.add(s);
    } else {
        throw except::invalid_argument("File::add: Type is not \"HEADER\" or \"FOOTER\"!");
    }
}

File File::copy() const {
    return File(*this);
}

File& File::operator=(const File& rhs) {
    protein_atoms = rhs.protein_atoms;
    hydration_atoms = rhs.hydration_atoms;
    header = rhs.header;
    footer = rhs.footer;
    terminate = rhs.terminate;
    return *this;
}

File& File::operator=(File&& rhs) = default;

void File::refresh() {
    if (protein_atoms.empty()) {return;}

    bool terminate_inserted = false;
    std::string chainID = "0"; int resSeq = 0; int serial = protein_atoms[0].serial;

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
        if (!terminate_inserted && a.get_type() == Record::RecordType::WATER) {
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

    chainID = protein_atoms[protein_atoms.size()-1].chainID;
    resSeq = protein_atoms[protein_atoms.size()-1].resSeq + 1;
    for (auto& a : hydration_atoms) {
        a.set_serial(serial++ % 100000);
        a.set_resSeq(resSeq++);
        a.set_chainID(chainID);
    }
}

std::unique_ptr<Reader> File::construct_reader(std::string path) {
    if (path.find(".xml") != std::string::npos || path.find(".XML") != std::string::npos) { // .xml file
        throw except::invalid_argument("File::construct_reader: .xml input files are not supported.");
    } else if (path.find(".pdb") != std::string::npos || path.find(".PDB") != std::string::npos) { // .pdb file
        return std::make_unique<PDBReader>(this);
    } else if (path.find(".ent") != std::string::npos || path.find(".ENT") != std::string::npos) { // .ent file
        return std::make_unique<PDBReader>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("File::construct_reader: Unsupported file extension of input file \"" + path + "\".");
    }
}

std::unique_ptr<Writer> File::construct_writer(std::string path) {
    if (path.find(".xml") != std::string::npos || path.find(".XML") != std::string::npos) { // .xml file
        throw except::invalid_argument("File::construct_writer: .xml input files are not supported.");
    } else if (path.find(".pdb") != std::string::npos || path.find(".PDB") != std::string::npos) { // .pdb file
        return std::make_unique<PDBWriter>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("File::construct_writer: Unsupported file extension of input file \"" + path + "\".");
    }
}