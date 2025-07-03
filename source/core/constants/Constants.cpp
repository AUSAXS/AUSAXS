// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <constants/Constants.h>
#include <constants/vdwTable.h>
#include <utility/StringUtils.h>
#include <io/ExistingFile.h>

#include <stdexcept>

using namespace ausaxs;

namespace ausaxs::constants {
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
    std::string symbols::carbon   = "C";
    std::string symbols::nitrogen = "N";
    std::string symbols::oxygen   = "O";

    const saxs::detail::SimpleMap<constants::atom_t> symbols::detail::string_to_atomt_map = std::unordered_map<std::string, constants::atom_t>{
        {"H", atom_t::H},   {"He", atom_t::He}, 
        {"Li", atom_t::Li}, {"Be", atom_t::Be}, {"B", atom_t::B},   {"C", atom_t::C},   {"N", atom_t::N}, {"O", atom_t::O}, {"F", atom_t::F},   {"Ne", atom_t::Ne}, 
        {"Na", atom_t::Na}, {"Mg", atom_t::Mg}, {"Al", atom_t::Al}, {"Si", atom_t::Si}, {"P", atom_t::P}, {"S", atom_t::S}, {"Cl", atom_t::Cl}, {"Ar", atom_t::Ar}, 
        {"K",  atom_t::K},  {"Ca", atom_t::Ca}, {"Sc", atom_t::Sc}, {"Ti", atom_t::Ti}, {"V", atom_t::V}, {"Cr", atom_t::Cr}, {"Mn", atom_t::Mn}, {"Fe", atom_t::Fe}, 
        {"Co", atom_t::Co}, {"Ni", atom_t::Ni}, {"Cu", atom_t::Cu}, {"Zn", atom_t::Zn}, {"Ga", atom_t::Ga}, {"Ge", atom_t::Ge}, {"As", atom_t::As}, {"Se", atom_t::Se}, 
        {"Br", atom_t::Br}, {"Kr", atom_t::Kr}, 
        {"Rb", atom_t::Rb}, {"Sr", atom_t::Sr}, {"Y",  atom_t::Y},  {"Zr", atom_t::Zr}, {"Nb", atom_t::Nb}, {"Mo", atom_t::Mo}, {"Tc", atom_t::Tc}, {"Ru", atom_t::Ru}, 
        {"Rh", atom_t::Rh}, {"Pd", atom_t::Pd}, {"Ag", atom_t::Ag}, {"Cd", atom_t::Cd}, {"In", atom_t::In}, {"Sn", atom_t::Sn}, {"Sb", atom_t::Sb}, {"Te", atom_t::Te},
        {"I",  atom_t::I},  {"Xe", atom_t::Xe}, 
        {"W",  atom_t::W},  {"M", atom_t::M}
    };

    residue::ResidueStorage hydrogen_atoms::residues;
}

constants::atomic_group_t constants::symbols::get_atomic_group(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom_type) {
    return hydrogen_atoms::residues.get_atomic_group(residue_name, atom_name, atom_type);
}

constants::atom_t constants::symbols::parse_element_string(const std::string& element_string) {
    return constants::symbols::detail::string_to_atomt_map.get(element_string);
}

constants::atomic_group_t constants::symbols::get_atomic_group(constants::atom_t atom_type, unsigned int hydrogens) {
    if (hydrogens == 0) {return constants::atomic_group_t::unknown;}
    switch (atom_type) {
        case constants::atom_t::C:
            switch (hydrogens) {
                case 1: return constants::atomic_group_t::CH;
                case 2: return constants::atomic_group_t::CH2;
                case 3: return constants::atomic_group_t::CH3;
                default: break;
            }
            break;
        case constants::atom_t::N:
            switch (hydrogens) {
                case 1: return constants::atomic_group_t::NH;
                case 2: return constants::atomic_group_t::NH2;
                case 3: return constants::atomic_group_t::NH3;
                default: break;
            }
            break;
        case constants::atom_t::O:
            switch (hydrogens) {
                case 1: return constants::atomic_group_t::OH;
                default: break;
            }
            break;
        case constants::atom_t::S:
            switch (hydrogens) {
                case 1: return constants::atomic_group_t::SH;
                default: break;
            }
            break;
        default: break;
    }
    return constants::atomic_group_t::unknown;
}

unsigned int constants::charge::nuclear::get_charge(atom_t atom) {
    switch(atom) {
        case atom_t::H:  return 1;
        case atom_t::He: return 2;

        case atom_t::Li: return 3;
        case atom_t::Be: return 4;
        case atom_t::B:  return 5;
        case atom_t::C:  return 6;
        case atom_t::N:  return 7;
        case atom_t::O:  return 8;
        case atom_t::F:  return 9;
        case atom_t::Ne: return 10;

        case atom_t::Na: return 11;
        case atom_t::Mg: return 12;
        case atom_t::Al: return 13;
        case atom_t::Si: return 14;
        case atom_t::P:  return 15;
        case atom_t::S:  return 16;
        case atom_t::Cl: return 17;
        case atom_t::Ar: return 18;

        case atom_t::K:  return 19;
        case atom_t::Ca: return 20;
        case atom_t::Sc: return 21;
        case atom_t::Ti: return 22;
        case atom_t::V:  return 23;
        case atom_t::Cr: return 24;
        case atom_t::Mn: return 25;
        case atom_t::Fe: return 26;
        case atom_t::Co: return 27;
        case atom_t::Ni: return 28;
        case atom_t::Cu: return 29;
        case atom_t::Zn: return 30;
        case atom_t::Ga: return 31;
        case atom_t::Ge: return 32;
        case atom_t::As: return 33;
        case atom_t::Se: return 34;
        case atom_t::Br: return 35;
        case atom_t::Kr: return 36;

        case atom_t::Rb: return 37;
        case atom_t::Sr: return 38;
        case atom_t::Y:  return 39;
        case atom_t::Zr: return 40;
        case atom_t::Nb: return 41;
        case atom_t::Mo: return 42;
        case atom_t::Tc: return 43;
        case atom_t::Ru: return 44;
        case atom_t::Rh: return 45;
        case atom_t::Pd: return 46;
        case atom_t::Ag: return 47;
        case atom_t::Cd: return 48;
        case atom_t::In: return 49;
        case atom_t::Sn: return 50;
        case atom_t::Sb: return 51;
        case atom_t::Te: return 52;
        case atom_t::I:  return 53;
        case atom_t::Xe: return 54;

        case atom_t::W:  return 74;
        case atom_t::M:  return 0;
        case atom_t::dummy: return 1;
        default: throw std::runtime_error("constants::charge::nuclear::get_charge: Unknown atom type \"" + constants::symbols::to_string(atom) + "\"");
    }
}

unsigned int constants::charge::ionic::get_charge(atom_t atom) {
    switch (atom) {
        case atom_t::Ca: return 2;
        case atom_t::Cl: return -1;
        case atom_t::Zn: return 2;
        default: throw std::runtime_error("constants::charge::ionic::get_charge: Unknown atom type \"" + constants::symbols::to_string(atom) + "\"");
    }
}

unsigned int constants::valence::get_valence(atom_t atom) {
    switch(atom) {
        case atom_t::H:  return 1;
        case atom_t::C:  return 4;
        case atom_t::N:  return 3;
        case atom_t::O:  return 2;
        case atom_t::F:  return 1;
        case atom_t::Ne: return 0;
        case atom_t::S:  return 2;
        case atom_t::P:  return 1;
        case atom_t::Cl: return 1;
        case atom_t::Fe: return 4;
        case atom_t::Se: return 2;
        case atom_t::M:  return 0;
        default: throw std::runtime_error("constants::valence::get_valence: Unknown atom type \"" + constants::symbols::to_string(atom) + "\"");
    }
}

double constants::radius::get_vdw_radius(atom_t atom) {
    switch(atom) {
        case atom_t::H:  return vdw::H;
        case atom_t::He: return vdw::He;
        case atom_t::Ne: return vdw::Ne;
        case atom_t::Ar: return vdw::Ar;

        case atom_t::Li: return vdw::Li;
        case atom_t::Be: return vdw::Be;
        case atom_t::B:  return vdw::B;
        case atom_t::C:  return vdw::C;
        case atom_t::N:  return vdw::N;
        case atom_t::O:  return vdw::O;
        case atom_t::F:  return vdw::F;

        case atom_t::Na: return vdw::Na;
        case atom_t::Mg: return vdw::Mg;
        case atom_t::Al: return vdw::Al;
        case atom_t::Si: return vdw::Si;
        case atom_t::P:  return vdw::P;
        case atom_t::S:  return vdw::S;
        case atom_t::Cl: return vdw::Cl;

        case atom_t::K:  return vdw::K;
        case atom_t::Ca: return vdw::Ca;
        case atom_t::Sc: return vdw::Sc;
        case atom_t::Ti: return vdw::Ti;
        case atom_t::V:  return vdw::V;
        case atom_t::Cr: return vdw::Cr;
        case atom_t::Mn: return vdw::Mn;
        case atom_t::Fe: return vdw::Fe;
        case atom_t::Co: return vdw::Co;
        case atom_t::Ni: return vdw::Ni;
        case atom_t::Cu: return vdw::Cu;
        case atom_t::Zn: return vdw::Zn;
        case atom_t::Ga: return vdw::Ga;
        case atom_t::Ge: return vdw::Ge;
        case atom_t::As: return vdw::As;
        case atom_t::Se: return vdw::Se;
        case atom_t::Br: return vdw::Br;

        case atom_t::Rb: return vdw::Rb;
        case atom_t::Sr: return vdw::Sr;
        case atom_t::Y:  return vdw::Y;
        case atom_t::Zr: return vdw::Zr;
        case atom_t::Nb: return vdw::Nb;
        case atom_t::Mo: return vdw::Mo;
        case atom_t::Tc: return vdw::Tc;
        case atom_t::Ru: return vdw::Ru;
        case atom_t::Rh: return vdw::Rh;
        case atom_t::Pd: return vdw::Pd;
        case atom_t::Ag: return vdw::Ag;
        case atom_t::Cd: return vdw::Cd;
        case atom_t::In: return vdw::In;
        case atom_t::Sn: return vdw::Sn;
        case atom_t::Sb: return vdw::Sb;
        case atom_t::Te: return vdw::Te;
        case atom_t::I:  return vdw::I;

        case atom_t::W:  return vdw::W;

        // fake elements
        case atom_t::M: return 0;
        case atom_t::dummy: {return radius::detail::dummy_radius;}
        default: throw std::runtime_error("constants::radius::get_vdw_radius: Unknown atom type \"" + constants::symbols::to_string(atom) + "\"");
    }
}

std::string ausaxs::constants::symbols::to_string(atom_t atom) {
    switch(atom) {
        case atom_t::H:  return "H";
        case atom_t::He: return "He";

        case atom_t::Li: return "Li";
        case atom_t::Be: return "Be";
        case atom_t::B:  return "B";
        case atom_t::C:  return "C";
        case atom_t::N:  return "N";
        case atom_t::O:  return "O";
        case atom_t::F:  return "F";
        case atom_t::Ne: return "Ne";

        case atom_t::Na: return "Na";
        case atom_t::Mg: return "Mg";
        case atom_t::Al: return "Al";
        case atom_t::Si: return "Si";
        case atom_t::P:  return "P";
        case atom_t::S:  return "S";
        case atom_t::Cl: return "Cl";
        case atom_t::Ar: return "Ar";

        case atom_t::K:  return "K";
        case atom_t::Ca: return "Ca";
        case atom_t::Sc: return "Sc";
        case atom_t::Ti: return "Ti";
        case atom_t::V:  return "V";
        case atom_t::Cr: return "Cr";
        case atom_t::Mn: return "Mn";
        case atom_t::Fe: return "Fe";
        case atom_t::Co: return "Co";
        case atom_t::Ni: return "Ni";
        case atom_t::Cu: return "Cu";
        case atom_t::Zn: return "Zn";
        case atom_t::Ga: return "Ga";
        case atom_t::Ge: return "Ge";
        case atom_t::As: return "As";
        case atom_t::Se: return "Se";
        case atom_t::Br: return "Br";

        case atom_t::Rb: return "Rb";
        case atom_t::Sr: return "Sr";
        case atom_t::Y:  return "Y";
        case atom_t::Zr: return "Zr";
        case atom_t::Nb: return "Nb";
        case atom_t::Mo: return "Mo";
        case atom_t::Tc: return "Tc";
        case atom_t::Ru: return "Ru";
        case atom_t::Rh: return "Rh";
        case atom_t::Pd: return "Pd";
        case atom_t::Ag: return "Ag";
        case atom_t::Cd: return "Cd";
        case atom_t::In: return "In";
        case atom_t::Sn: return "Sn";
        case atom_t::Sb: return "Sb";
        case atom_t::Te: return "Te";
        case atom_t::I:  return "I";
        case atom_t::Xe: return "Xe";

        case atom_t::W:  return "W";
        case atom_t::M:  return "M";
        case atom_t::dummy: return "#";
        default: throw std::runtime_error("constants::symbols::to_string: Unknown atom type \"" + std::to_string(static_cast<int>(atom)) + "\"");
    }
}

std::string ausaxs::constants::symbols::to_string(atomic_group_t group) {
    switch(group) {
        case atomic_group_t::CH: return "CH";
        case atomic_group_t::CH2: return "CH2";
        case atomic_group_t::CH3: return "CH3";
        case atomic_group_t::NH: return "NH";
        case atomic_group_t::NH2: return "NH2";
        case atomic_group_t::NH3: return "NH3";
        case atomic_group_t::OH: return "OH";
        case atomic_group_t::SH: return "SH";
        case atomic_group_t::unknown: return "unknown";
        default: throw std::runtime_error("constants::symbols::to_string: Unknown atomic group \"" + std::to_string(static_cast<int>(group)) + "\"");
    }
}