#include <constants/Constants.h>
#include <utility/StringUtils.h>
#include <io/ExistingFile.h>

namespace constants {
    filetypes::detail::FileType::FileType(std::vector<std::string> extensions) {
        for (auto& ext : extensions) {
            ext = utility::to_lowercase(ext);
        }
        this->extensions = extensions;
    }

    bool filetypes::detail::FileType::validate(const io::ExistingFile& path) const {
        std::string file_ext = utility::to_lowercase(path.extension()); 
        for (auto ext : extensions) {
            if (file_ext == ext) {
                return true;
            }
        }
        return false;
    }

    const saxs::detail::SimpleMap<char> name_1symbol_map = std::unordered_map<std::string, char>{
        {"glycine", 'G'}, {"alanine", 'A'}, {"valine", 'V'}, {"leucine", 'L'}, {"isoleucine", 'I'}, {"phenylalanine", 'F'}, {"tyrosine", 'Y'}, 
        {"tryptophan", 'W'}, {"aspartic_acid", 'D'}, {"glutamic_acid", 'E'}, {"serine", 'S'}, {"threonine", 'T'}, {"asparagine", 'N'}, 
        {"glutamine", 'Q'}, {"lysine", 'K'}, {"arginine", 'R'}, {"histidine", 'H'}, {"methionine", 'M'}, {"cysteine", 'C'}, {"proline", 'P'}
    };

    const saxs::detail::SimpleMap<std::string> name_3symbol_map = std::unordered_map<std::string, std::string>{
        {"glycine", "GLY"}, {"alanine", "ALA"}, {"valine", "VAL"}, {"leucine", "LEU"}, {"isoleucine", "ILE"}, {"phenylalanine", "PHE"}, {"tyrosine", "TYR"}, 
        {"tryptophan", "TRP"}, {"aspartic_acid", "ASP"}, {"glutamic_acid", "GLU"}, {"serine", "SER"}, {"threonine", "THR"}, {"asparagine", "ASN"}, 
        {"glutamine", "GLN"}, {"lysine", "LYS"}, {"arginine", "ARG"}, {"histidine", "HIS"}, {"methionine", "MET"}, {"cysteine", "CYS"}, {"proline", "PRO"}
    };

    const saxs::detail::SimpleMap<double> volume::amino_acids = std::unordered_map<std::string, double>{
        {"GLY", 66.4}, {"ALA", 91.5}, {"VAL", 141.7}, {"LEU", 167.9}, {"ILE", 168.8}, {"PHE", 203.5}, {"TYR", 203.6}, {"TRP", 237.6}, 
        {"ASP", 113.6}, {"GLU", 140.6}, {"SER", 99.1}, {"THR", 122.1}, {"ASN", 135.2}, {"GLN", 161.1}, {"LYS", 176.2}, {"ARG", 180.8}, 
        {"HIS", 167.3}, {"MET", 170.8}, {"CYS", 105.6}, {"PRO", 129.3}
    };

    double radius::detail::dummy_radius = 1;
    void radius::set_dummy_radius(double radius) {
        radius::detail::dummy_radius = radius;
    }

    std::string symbols::hydrogen = "H";
    std::string symbols::carbon = "C";
    std::string symbols::nitrogen = "N";
    std::string symbols::oxygen = "O";

    const saxs::detail::SimpleMap<constants::atom_t> symbols::detail::string_to_atomt_map = std::unordered_map<std::string, constants::atom_t>{
        {"e", atom_t::e}, {"H", atom_t::H}, {"He", atom_t::He}, {"Li", atom_t::Li}, {"Be", atom_t::Be}, {"B", atom_t::B}, {"C", atom_t::C}, {"N", atom_t::N},
        {"O", atom_t::O}, {"F", atom_t::F}, {"Ne", atom_t::Ne}, {"Na", atom_t::Na}, {"Mg", atom_t::Mg}, {"Al", atom_t::Al}, {"Si", atom_t::Si}, {"P", atom_t::P},
        {"S", atom_t::S}, {"Cl", atom_t::Cl}, {"Ar", atom_t::Ar}, {"K", atom_t::K}, {"Ca", atom_t::Ca}, {"Sc", atom_t::Sc}, {"Ti", atom_t::Ti}, {"V", atom_t::V},
        {"Cr", atom_t::Cr}, {"Mn", atom_t::Mn}, {"Fe", atom_t::Fe}, {"Co", atom_t::Co}, {"Ni", atom_t::Ni}, {"Cu", atom_t::Cu}, {"Zn", atom_t::Zn}, {"W", atom_t::W},
        {"M", atom_t::M}
    };

    //* note: this must be initialized *after* symbols::detail::string_to_atomt_map
    parser::residue::ResidueStorage hydrogen_atoms::residues;
}

constants::atom_t constants::symbols::parse_element_string(const std::string& element_string) {
    return constants::symbols::detail::string_to_atomt_map.get(element_string);
}

std::string constants::symbols::write_element_string(atom_t atom) {
    switch(atom) {
        case atom_t::e: return "e";
        case atom_t::H: return "H";
        case atom_t::He: return "He";
        case atom_t::Li: return "Li";
        case atom_t::Be: return "Be";
        case atom_t::B: return "B";
        case atom_t::C: return "C";
        case atom_t::N: return "N";
        case atom_t::O: return "O";
        case atom_t::F: return "F";
        case atom_t::Ne: return "Ne";
        case atom_t::Na: return "Na";
        case atom_t::Mg: return "Mg";
        case atom_t::Al: return "Al";
        case atom_t::Si: return "Si";
        case atom_t::P: return "P";
        case atom_t::S: return "S";
        case atom_t::Cl: return "Cl";
        case atom_t::Ar: return "Ar";
        case atom_t::K: return "K";
        case atom_t::Ca: return "Ca";
        case atom_t::Sc: return "Sc";
        case atom_t::Ti: return "Ti";
        case atom_t::V: return "V";
        case atom_t::Cr: return "Cr";
        case atom_t::Mn: return "Mn";
        case atom_t::Fe: return "Fe";
        case atom_t::Co: return "Co";
        case atom_t::Ni: return "Ni";
        case atom_t::Cu: return "Cu";
        case atom_t::Zn: return "Zn";
        case atom_t::W: return "W";
        case atom_t::M: return "M";
        case atom_t::dummy: return "dummy";
        default: throw std::runtime_error("constants::symbols::write_element_string: Unknown atom type");
    }
}

unsigned int constants::charge::get_charge(atom_t atom) {
    switch(atom) {
        case atom_t::e: return -1;
        case atom_t::H: return 1;
        case atom_t::He: return 2;
        case atom_t::Li: return 3;
        case atom_t::Be: return 4;
        case atom_t::B: return 5;
        case atom_t::C: return 6;
        case atom_t::N: return 7;
        case atom_t::O: return 8;
        case atom_t::F: return 9;
        case atom_t::Ne: return 10;
        case atom_t::Na: return 11;
        case atom_t::Mg: return 12;
        case atom_t::Al: return 13;
        case atom_t::Si: return 14;
        case atom_t::P: return 15;
        case atom_t::S: return 16;
        case atom_t::Cl: return 17;
        case atom_t::Ar: return 18;
        case atom_t::K: return 19;
        case atom_t::Ca: return 20;
        case atom_t::Sc: return 21;
        case atom_t::Ti: return 22;
        case atom_t::V: return 23;
        case atom_t::Cr: return 24;
        case atom_t::Mn: return 25;
        case atom_t::Fe: return 26;
        case atom_t::Co: return 27;
        case atom_t::Ni: return 28;
        case atom_t::Cu: return 29;
        case atom_t::Zn: return 30;
        case atom_t::W: return 74;
        case atom_t::M: return 0;
        case atom_t::dummy: return 0;
        default: throw std::runtime_error("constants::charge::get_charge: Unknown atom type");
    }
}

unsigned int constants::valence::get_valence(atom_t atom) {
    switch(atom) {
        case atom_t::H: return 1;
        case atom_t::C: return 4;
        case atom_t::N: return 3;
        case atom_t::O: return 2;
        case atom_t::F: return 1;
        case atom_t::Ne: return 0;
        case atom_t::S: return 2;
        case atom_t::P: return 1;
        case atom_t::Cl: return 1;
        case atom_t::M: return 0;
        default: throw std::runtime_error("constants::valence::get_valence: Unknown atom type");
    }
}

double constants::mass::get_mass(atom_t atom) {
    switch(atom) {
        case atom_t::H: return 1.0079;
        case atom_t::He: return 4.0026;
        case atom_t::Li: return 6.941;
        case atom_t::Be: return 9.0122;
        case atom_t::B: return 10.811;
        case atom_t::C: return 12.0107;
        case atom_t::N: return 14.0067;
        case atom_t::O: return 15.9994;
        case atom_t::F: return 18.9984;
        case atom_t::Ne: return 20.1797;
        case atom_t::Na: return 22.9897;
        case atom_t::Mg: return 24.305;
        case atom_t::Al: return 26.9815;
        case atom_t::Si: return 28.0855;
        case atom_t::P: return 30.9738;
        case atom_t::S: return 32.065;
        case atom_t::Cl: return 35.453;
        case atom_t::Ar: return 39.948;
        case atom_t::K: return 39.0983;
        case atom_t::Ca: return 40.078;
        case atom_t::Sc: return 44.9559;
        case atom_t::Ti: return 47.867;
        case atom_t::V: return 50.9415;
        case atom_t::Cr: return 51.9961;
        case atom_t::Mn: return 54.938;
        case atom_t::Fe: return 55.845;
        case atom_t::Co: return 58.9332;
        case atom_t::Ni: return 58.6934;
        case atom_t::Cu: return 63.546;
        case atom_t::Zn: return 65.39;
        case atom_t::W: return 183.84;
        case atom_t::M: return 0;
        case atom_t::dummy: return 1;
        default: throw std::runtime_error("constants::mass::get_mass: Unknown atom type");   
    }
}

double constants::radius::get_vdw_radius(atom_t atom) {
    // switch(atom) {
    //     case atom_t::O: return 1.5;
    //     default: return 2.4;
    // }
    switch(atom) {
        // Rowland 1996, crystallographic
        case atom_t::H: return 1.1;

        // wikipedia, likely also crystallographic
        case atom_t::He: return 1.4;
        case atom_t::Ne: return 1.54;
        case atom_t::Ar: return 1.88;

        // Batsanov, 2001, equilibrium values
        case atom_t::Li: return 2.63;
        case atom_t::Be: return 2.23;
        case atom_t::B: return 2.05;
        case atom_t::C: return 1.96;
        case atom_t::N: return 1.79;
        case atom_t::O: return 1.71;
        case atom_t::F: return 1.65;
        case atom_t::Na: return 2.77;
        case atom_t::Mg: return 2.42;
        case atom_t::Al: return 2.40;
        case atom_t::Si: return 2.26;
        case atom_t::P: return 2.14;
        case atom_t::S: return 2.06;
        case atom_t::Cl: return 2.05;
        case atom_t::K: return 3.02;
        case atom_t::Ca: return 2.78;
        case atom_t::Sc: return 2.62;
        case atom_t::Ti: return 2.44;
        case atom_t::V: return 2.27;
        case atom_t::Cr: return 2.23;
        case atom_t::Mn: return 2.25;
        case atom_t::Fe: return 2.27;
        case atom_t::Co: return 2.25;
        case atom_t::Ni: return 2.23;
        case atom_t::Cu: return 2.27;
        case atom_t::Zn: return 2.24;
        case atom_t::W: return 2.36;

        // fake elements
        case atom_t::M: return 0;
        case atom_t::dummy: {
            std::cout << "constants::radius::get_vdw_radius: Warning: dummy atom radius requested. Returning " << radius::detail::dummy_radius << " Ångström." << std::endl;
            return radius::detail::dummy_radius;
        }

        default: throw std::runtime_error("constants::radius::get_vdw_radius: Unknown atom type");
    }
}