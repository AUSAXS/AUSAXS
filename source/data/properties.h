#pragma once

#include <map>
#include <string>

using std::string;

namespace property {
    const std::map<string, char> name_1symbol_map = {{"glycine", 'G'}, {"alanine", 'A'}, {"valine", 'V'}, {"leucine", 'L'}, {"isoleucine", 'I'}, 
        {"phenylalanine", 'F'}, {"tyrosine", 'Y'}, {"tryptophan", 'W'}, {"aspartic_acid", 'D'}, {"glutamic_acid", 'E'}, {"serine", 'S'}, 
        {"threonine", 'T'}, {"asparagine", 'N'}, {"glutamine", 'Q'}, {"lysine", 'K'}, {"arginine", 'R'}, {"histidine", 'H'}, {"methionine", 'M'}, 
        {"cysteine", 'C'}, {"proline", 'P'}
    };

    const std::map<string, string> name_3symbol_map = {{"glycine", "GLY"}, {"alanine", "ALA"}, {"valine", "VAL"}, {"leucine", "LEU"}, 
        {"isoleucine", "ILE"}, {"phenylalanine", "PHE"}, {"tyrosine", "TYR"}, {"tryptophan", "TRP"}, {"aspartic_acid", "ASP"}, {"glutamic_acid", "GLU"}, 
        {"serine", "SER"}, {"threonine", "THR"}, {"asparagine", "ASN"}, {"glutamine", "GLN"}, {"lysine", "LYS"}, {"arginine", "ARG"}, {"histidine", "HIS"}, 
        {"methionine", "MET"}, {"cysteine", "CYS"}, {"proline", "PRO"}
    };


    // taken from https://doi.org/10.1088/0034-4885/39/10/001 
    // all values are in Å^3
    namespace volume {
        constexpr double glycine = 66.4;
        constexpr double alanine = 91.5;
        constexpr double valine = 141.7;
        constexpr double leucine = 167.9;
        constexpr double isoleucine = 168.8;
        constexpr double phenylalanine = 203.5;
        constexpr double tyrosine = 203.6;
        constexpr double tryptophan = 237.6;
        constexpr double aspartic_acid = 113.6;
        constexpr double glutamic_acid = 140.6;
        constexpr double serine = 99.1;
        constexpr double threonine = 122.1;
        constexpr double asparagine = 135.2;
        constexpr double glutamine = 161.1;
        constexpr double lysine = 176.2;
        constexpr double arginine = 180.8;
        constexpr double histidine = 167.3;
        constexpr double methionine = 170.8;
        constexpr double cysteine = 105.6;
        constexpr double proline = 129.3;

        // get the volume of a 3symbol amino acid
        const std::map<string, double> get = {{"GLY", glycine}, {"ALA", alanine}, {"VAL", valine}, {"LEU", leucine}, {"ILE", isoleucine}, 
            {"PHE", phenylalanine}, {"TYR", tyrosine}, {"TRP", tryptophan}, {"ASP", aspartic_acid}, {"GLU", glutamic_acid}, {"SER", serine}, 
            {"THR", threonine}, {"ASN", asparagine}, {"GLN", glutamine}, {"LYS", lysine}, {"ARG", arginine}, {"HIS", histidine}, 
            {"MET", methionine}, {"CYS", cysteine}, {"PRO", proline},
        };
    }

    // atomic weights taken from https://www.britannica.com/science/atomic-weight
    namespace weight {
        constexpr double H = 1.01;
        constexpr double He = 4.00;
        constexpr double Li = 6.95;
        constexpr double C = 12.01;
        constexpr double N = 14.01;
        constexpr double O = 16.00;
        constexpr double S = 32.06;

        // get the weight of an atom
        const std::map<string, double> atomic = {{"H", H}, {"He", He}, {"Li", Li}, {"C", C}, {"N", N}, {"O", O}, {"S", S}};
    }

    namespace charge {
        constexpr int H = 1;
        constexpr int He = 2;
        constexpr int Li = 3;
        constexpr int C = 6;
        constexpr int N = 7;
        constexpr int O = 8;
        constexpr int S = 16;

        // get the charge Z of an atom
        const std::map<string, int> get = {{"H", H}, {"He", He}, {"Li", Li}, {"C", C}, {"N", N}, {"O", O}, {"S", S}};
    }
}