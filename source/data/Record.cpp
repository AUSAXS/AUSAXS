#include <data/Record.h>

#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

RecordType Record::get_type(const std::string& s) {
    auto str = utility::remove_all(s, " \r"); // remove any space or carriage returns, since programs are inconsistent with the spacing after e.g. END or TER
    if (type_map.count(str) == 1) {
        return type_map.at(str);
    }
    throw except::parse_error("Record::get_type: Could not determine type \"" + s + "\"");
}

// Maps PDB types to a Record. Effectively determines how they are treated by the code.
//      ATOMs will be converted to an Atom object. 
//      TERMINATE will be parsed as a Terminate statement. 
//      All HEADERs will be combined into a single string. HEADERs must always be above the ATOM section.
//      FOOTERs are treated identically to HEADERs, but must always be after the ATOM section.
//      NOTYPEs are ignored. 
const std::unordered_map<std::string, RecordType> Record::type_map = {
    {"ATOM"  , RecordType::ATOM}, {"HETATM", RecordType::ATOM},

    {"TER"   , RecordType::TERMINATE}, 

    {"HEADER", RecordType::HEADER}, {"TITLE" , RecordType::HEADER}, {"COMPND", RecordType::HEADER}, {"SOURCE", RecordType::HEADER}, 
    {"KEYWDS", RecordType::HEADER}, {"EXPDTA", RecordType::HEADER}, {"AUTHOR", RecordType::HEADER}, {"REVDAT", RecordType::HEADER}, 
    {"JRNL"  , RecordType::HEADER}, {"REMARK", RecordType::HEADER}, {"DBREF" , RecordType::HEADER}, {"SEQRES", RecordType::HEADER}, 
    {"FORMUL", RecordType::HEADER}, {"HELIX" , RecordType::HEADER}, {"SHEET" , RecordType::HEADER}, {"SSBOND", RecordType::HEADER}, 
    {"CRYST1", RecordType::HEADER}, {"ORIGX1", RecordType::HEADER}, {"ORIGX2", RecordType::HEADER}, {"ORIGX3", RecordType::HEADER}, 
    {"SCALE1", RecordType::HEADER}, {"SCALE2", RecordType::HEADER}, {"SCALE3", RecordType::HEADER}, {"HET"   , RecordType::HEADER}, 
    {"HETNAM", RecordType::HEADER}, {"HETSYN", RecordType::HEADER}, {"FORMUL", RecordType::HEADER}, {"CISPEP", RecordType::HEADER}, 
    {"SITE"  , RecordType::HEADER}, {"SEQADV", RecordType::HEADER}, {"LINK"  , RecordType::HEADER}, {"MODEL" , RecordType::HEADER}, 
    {"LINKR" , RecordType::HEADER}, {"SPRSDE", RecordType::HEADER}, {"MODRES", RecordType::HEADER}, {"MTRIX1", RecordType::HEADER}, 
    {"MTRIX2", RecordType::HEADER}, {"MTRIX3", RecordType::HEADER}, {"DBREF1", RecordType::HEADER}, {"DBREF2", RecordType::HEADER},
    {"CAVEAT", RecordType::HEADER},

    {"CONECT", RecordType::FOOTER}, {"MASTER", RecordType::FOOTER}, {"END"   , RecordType::FOOTER}, {"ENDMDL", RecordType::FOOTER},

    {"ANISOU", RecordType::NOTYPE}
};