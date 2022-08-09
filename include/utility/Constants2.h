#pragma once

#include <map>
#include <string>

constexpr double simple_pow(double val, unsigned int power) {
    double sum = 1;
    for (unsigned int i = 0; i < power; i++) {
        sum *= val;
    }
    return sum;
}

template<typename K, typename V>
struct Storage {
    /**
     * @brief Create a new empty storage. 
     */
    Storage() {}

    /**
     * @brief Create a new storage from a map.
     */
    Storage(std::map<K, V> map) : data(map) {}

    /**
     * @brief Get a value from the storage. 
     */
    V get(K key) {
        if (data.find(key) == data.end()) {
            throw except::map_error("Key " + key + " not found in map");
        }
        return data.at(key);
    }

    /**
     * @brief Insert a key-value pair into the storage.
     */
    void insert(K key, V val) {
        data.emplace(key, val);
    }

    /**
     * @brief Get a value from the storage. 
     *        The key is converted to lowercase. 
     */
    template<>
    V get<std::string>(std::string key) {
        std::string k2 = utility::to_lowercase(key);
        if (data.find(k2) == data.end()) {
            throw except::map_error("Key " + k2 + " not found in map");
        }
        return data.at(k2);        
    }

    /**
     * @brief Insert a key-value pair into the storage. 
     *        The key is converted to lowercase. 
     */
    template<>
    void insert<std::string> (std::string key, V val) {
        std::string k2 = utility::to_lowercase(key);
        data.emplace(k2, val);
    }

    std::map<K, V> data;
};

/**
 * @brief \namespace constants
 * 
 * This namespace contains all constants used in this project. 
 */
namespace constants {
    /**
     * @brief \namespace radius
     * 
     * This namespace contains all the radius constants used in this project. 
     */
    namespace radius {
        constexpr double electron = 0.0000281794; // electron radius in units of Ångström
    }
    constexpr double Avogadro = 6.02214076e-23; // mol^-1

    /**
     * @brief \namespace Relative units.
     * 
     * This namespace contains all the unit conversion constants used in this project. 
     * The basic units are for:
     *     mass: Dalton
     *     length: Å
     *     charge: e
     */
    namespace unit {
        constexpr double gm = 1.66054e-24; // Dalton --> grams
        constexpr double mg = 1.66054e-21; // Dalton --> mg
        constexpr double cm = 1e-8; // Ångström --> cm

        constexpr double mL = simple_pow(unit::cm, 3); // Ångström^3 --> mL
    }

    /**
     * @brief \namespace Absolute units.
     * 
     * This namespace contains all the absolute unit conversion constants. 
     */
    namespace SI {
        namespace mass {
            constexpr double gm = 1e-3;
            constexpr double mg = 1e-6;
            constexpr double u = 1.66053*1e-27;
        }

        namespace length {
            constexpr double cm = 1e-3;
            constexpr double A = 1e-8; // Ångström
        }

        namespace volume {
            constexpr double A3 = 1e-24; // Ångström^3
            constexpr double cm3 = 1e-9;
        }
    }

    // The 1-symbol names of all amino acids. 
    const std::map<std::string, char> name_1symbol_map = {{"glycine", 'G'}, {"alanine", 'A'}, {"valine", 'V'}, {"leucine", 'L'}, {"isoleucine", 'I'}, 
        {"phenylalanine", 'F'}, {"tyrosine", 'Y'}, {"tryptophan", 'W'}, {"aspartic_acid", 'D'}, {"glutamic_acid", 'E'}, {"serine", 'S'}, 
        {"threonine", 'T'}, {"asparagine", 'N'}, {"glutamine", 'Q'}, {"lysine", 'K'}, {"arginine", 'R'}, {"histidine", 'H'}, {"methionine", 'M'}, 
        {"cysteine", 'C'}, {"proline", 'P'}
    };
    const Storage name_1symbol_map = std::map<std::string, char>{{
        {"glycine", 'G'}, {"alanine", 'A'}, {"valine", 'V'}, {"leucine", 'L'}, {"isoleucine", 'I'}, {"phenylalanine", 'F'}, {"tyrosine", 'Y'}, 
        {"tryptophan", 'W'}, {"aspartic_acid", 'D'}, {"glutamic_acid", 'E'}, {"serine", 'S'}, {"threonine", 'T'}, {"asparagine", 'N'}, 
        {"glutamine", 'Q'}, {"lysine", 'K'}, {"arginine", 'R'}, {"histidine", 'H'}, {"methionine", 'M'}, {"cysteine", 'C'}, {"proline", 'P'}
    }};

    // The 3-symbol names of all amino acids. 
    const Storage name_3symbol_map = std::map<std::string, std::string>{{
        {"glycine", "GLY"}, {"alanine", "ALA"}, {"valine", "VAL"}, {"leucine", "LEU"}, {"isoleucine", "ILE"}, {"phenylalanine", "PHE"}, {"tyrosine", "TYR"}, 
        {"tryptophan", "TRP"}, {"aspartic_acid", "ASP"}, {"glutamic_acid", "GLU"}, {"serine", "SER"}, {"threonine", "THR"}, {"asparagine", "ASN"}, 
        {"glutamine", "GLN"}, {"lysine", "LYS"}, {"arginine", "ARG"}, {"histidine", "HIS"}, {"methionine", "MET"}, {"cysteine", "CYS"}, {"proline", "PRO"}
    }};


    /**
     * @brief \namespace volume
     * 
     * This namespace contains the volume of all amino acids. 
     * They are taken from https://doi.org/10.1088/0034-4885/39/10/001.
     * All values are in Å^3
     */
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
        const Storage get = std::map<std::string, double>{{
            {"GLY", glycine}, {"ALA", alanine}, {"VAL", valine}, {"LEU", leucine}, {"ILE", isoleucine}, {"PHE", phenylalanine}, {"TYR", tyrosine}, 
            {"TRP", tryptophan}, {"ASP", aspartic_acid}, {"GLU", glutamic_acid}, {"SER", serine}, {"THR", threonine}, {"ASN", asparagine}, 
            {"GLN", glutamine}, {"LYS", lysine}, {"ARG", arginine}, {"HIS", histidine}, {"MET", methionine}, {"CYS", cysteine}, {"PRO", proline}
        }};
    }

    /**
     * @brief \namespace mass 
     * 
     * This namespace contains the masses of the most common atomic elements encountered in SAXS. 
     */
    namespace mass {
        constexpr double H = 1.0079;
        constexpr double He = 4.0026;
        constexpr double Li = 6.941;
        constexpr double Be = 9.0122;
        constexpr double B = 10.811;
        constexpr double C = 12.0107;
        constexpr double N = 14.0067;
        constexpr double O = 15.9994;
        constexpr double F = 18.9984;
        constexpr double Ne = 20.1797;
        constexpr double Na = 22.9897;
        constexpr double Mg = 24.305;
        constexpr double Al = 26.9815;
        constexpr double Si = 28.0855;
        constexpr double P = 30.9738;
        constexpr double S = 32.065;
        constexpr double Cl = 35.453;
        constexpr double Ar = 39.948;
        constexpr double K = 39.0983;
        constexpr double Ca = 40.078;
        constexpr double Sc = 44.9559;
        constexpr double Ti = 47.867;
        constexpr double V = 50.9415;
        constexpr double Cr = 51.9961;
        constexpr double Mn = 54.938;
        constexpr double Fe = 55.845;
        constexpr double Co = 58.9332;
        constexpr double Ni = 58.6934;
        constexpr double Cu = 63.546;
        constexpr double Zn = 65.39;
        constexpr double W = 183.84;
        
        // get the weight of an atom
        const std::map<std::string, double> atomic = {{"H", H}, {"He", He}, {"HE", He}, {"Li", Li}, {"LI", Li}, {"Be", Be}, {"BE", Be}, {"B", B}, {"C", C}, 
            {"N", N}, {"O", O}, {"F", F}, {"Ne", Ne}, {"NE", Ne}, {"Na", Na}, {"NA", Na}, {"Mg", Mg}, {"MG", Mg}, {"Al", Al}, {"AL", Al}, {"Si", Si}, 
            {"SI", Si}, {"P", P}, {"S", S}, {"Cl", Cl}, {"CL", Cl}, {"Ar", Ar}, {"AR", Ar}, {"K", K}, {"Ca", Ca}, {"CA", Ca}, {"Sc", Sc}, {"SC", Sc}, 
            {"Ti", Ti}, {"TI", Ti}, {"V", V}, {"Cr", Cr}, {"CR", Cr}, {"Mn", Mn}, {"MN", Mn}, {"Fe", Fe}, {"FE", Fe}, {"Co", Co}, {"CO", Co}, {"Ni", Ni}, 
            {"NI", Ni}, {"Cu", Cu}, {"CU", Cu}, {"Zn", Zn}, {"ZN", Zn}, {"W", W}
        };
        namespace density {
            constexpr double water = 0.9982067*SI::mass::u/SI::volume::A3; // u/Å^3
        }
    }

    /**
     * @brief \namespace charge
     * 
     * This namespace contains the net charge of the most common atomic elements encountered in SAXS. 
     */
    namespace charge {
        constexpr unsigned int e = 1;
        constexpr unsigned int H = 1;
        constexpr unsigned int He = 2;
        constexpr unsigned int Li = 3;
        constexpr unsigned int Be = 4;
        constexpr unsigned int B = 5;
        constexpr unsigned int C = 6;
        constexpr unsigned int N = 7;
        constexpr unsigned int O = 8;
        constexpr unsigned int F = 9;
        constexpr unsigned int Ne = 10;
        constexpr unsigned int Na = 11;
        constexpr unsigned int Mg = 12;
        constexpr unsigned int Al = 13;
        constexpr unsigned int Si = 14;
        constexpr unsigned int P = 15;
        constexpr unsigned int S = 16;
        constexpr unsigned int Cl = 17;
        constexpr unsigned int Ar = 18;
        constexpr unsigned int K = 19;
        constexpr unsigned int Ca = 20;
        constexpr unsigned int Sc = 21;
        constexpr unsigned int Ti = 22;
        constexpr unsigned int V = 23;
        constexpr unsigned int Cr = 24;
        constexpr unsigned int Mn = 25;
        constexpr unsigned int Fe = 26;
        constexpr unsigned int Co = 27;
        constexpr unsigned int Ni = 28;
        constexpr unsigned int Cu = 29;
        constexpr unsigned int Zn = 30;
        constexpr unsigned int W = 74;

        // get the charge Z of an atom
        const std::map<std::string, unsigned int> atomic = {{"e", e}, {"H", H}, {"He", He}, {"HE", He}, {"Be", Be}, {"BE", Be}, {"B", B}, {"C", C}, 
            {"N", N}, {"O", O}, {"F", F}, {"Ne", Ne}, {"NE", Ne}, {"Na", Na}, {"NA", Na}, {"Mg", Mg}, {"MG", Mg}, {"Al", Al}, {"AL", Al}, {"Si", Si}, 
            {"SI", Si}, {"P", P}, {"S", S}, {"Cl", Cl}, {"CL", Cl}, {"Ar", Ar}, {"AR", Ar}, {"K", K}, {"Ca", Ca}, {"CA", Ca}, {"Sc", Sc}, {"SC", Sc}, 
            {"Ti", Ti}, {"TI", Ti}, {"V", V}, {"Cr", Cr}, {"CR", Cr}, {"Mn", Mn}, {"MN", Mn}, {"Fe", Fe}, {"FE", Fe}, {"Co", Co}, {"CO", Co}, {"Ni", Ni}, 
            {"NI", Ni}, {"Cu", Cu}, {"CU", Cu}, {"Zn", Zn}, {"ZN", Zn}, {"W", W}
        };

        namespace density {
            constexpr double water = 0.334; // e/Å^3
        }
    }

    /**
     * @brief valence
     */
    namespace valence {
        constexpr unsigned int H = 1;
        constexpr unsigned int C = 4;
        constexpr unsigned int N = 3;
        constexpr unsigned int O = 2;
        constexpr unsigned int F = 1;
        constexpr unsigned int Ne = 0;
        constexpr unsigned int S = 2;
        constexpr unsigned int P = 1;
        constexpr unsigned int Cl = 1;

        // get the valence of an atom
        const std::map<std::string, unsigned int> get = {{"H", H}, {"C", C}, {"N", N}, {"O", O}, {"F", F}, {"NE", Ne}, {"Ne", Ne}, {"S", S}, {"P", P}, 
        {"CL", Cl}, {"Cl", Cl}};
    }

    /**
     * @brief \namespace hydrogen_atoms
     * 
     * This namespace contains the number of hydrogen atoms attached to all amino acids. 
     */
    namespace hydrogen_atoms {
        namespace none {
            const std::map<std::string, int> get = {{"", 0}};
        }
        namespace water {
            constexpr int O = 2;
            const std::map<std::string, int> get = {{"O", O}};
        }
        namespace glycine {
            constexpr int N = 1;
            constexpr int CA = 2;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}};
        }
        namespace alanine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}};
        }
        namespace valine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 1;
            constexpr int CG1 = 3;
            constexpr int CG2 = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG1", CG1}, {"CG2", CG2}};
        }
        namespace leucine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 1;
            constexpr int CD1 = 3;
            constexpr int CD2 = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, {"CD2", CD2}};
        }
        namespace isoleucine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 1;
            constexpr int CG2 = 3;
            constexpr int CG1 = 2;
            constexpr int CD1 = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG2", CG2}, {"CG1", CG1}, {"CD1", CD1}};
        }
        namespace phenylalanine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int CD1 = 1;
            constexpr int CD2 = 1;
            constexpr int CE1 = 1;
            constexpr int CE2 = 1;
            constexpr int CZ = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
                {"CD2", CD2}, {"CE1", CD1}, {"CE2", CD2}, {"CZ", CZ}};
        }
        namespace tyrosine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int CD1 = 1;
            constexpr int CD2 = 1;
            constexpr int CE1 = 1;
            constexpr int CE2 = 1;
            constexpr int CZ = 0;
            constexpr int OH = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
                {"CD2", CD2}, {"CE1", CE1}, {"CE2", CE2}, {"CZ", CZ}, {"OH", OH}};
        }
        namespace tryptophan {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int CD1 = 1;
            constexpr int CD2 = 0;
            constexpr int NE1 = 1;
            constexpr int CE2 = 0;
            constexpr int CE3 = 1;
            constexpr int CZ2 = 1;
            constexpr int CZ3 = 1;
            constexpr int CH2 = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
                {"CD2", CD2}, {"NE1", NE1}, {"CE2", CE2}, {"CE3", CE3}, {"CZ2", CZ2}, {"CZ3", CZ3}, {"CH2", CH2}};
        }
        namespace aspartic_acid {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int OD1 = 0;
            constexpr int OD2 = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"OD1", OD1}, 
                {"OD2", OD2}};
        }
        namespace glutamic_acid {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int CD = 0;
            constexpr int OE1 = 0;
            constexpr int OE2 = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
                {"OE1", OE1}, {"OE2", OE2}};
        }
        namespace serine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int OG = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"OG", OG}};
        }
        namespace threonine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 1;
            constexpr int OG1 = 1;
            constexpr int CG2 = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"OG1", OG1}, {"CG2", CG2}};
        }
        namespace asparagine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int OD1 = 0;
            constexpr int ND2 = 2;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"OD1", OD1}, 
                {"ND2", ND2}};
        }
        namespace glutamine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int CD = 0;
            constexpr int OE1 = 0;
            constexpr int NE2 = 2;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
                {"OE1", OE1}, {"NE2", NE2}};
        }
        namespace lysine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int CD = 2;
            constexpr int CE = 2;
            constexpr int NZ = 2;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
                {"CE", CE}, {"NZ", NZ}};
        }
        namespace arginine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int CD = 2;
            constexpr int NE = 1;
            constexpr int CZ = 0;
            constexpr int NH1 = 1;
            constexpr int NH2 = 2;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
                {"NE", NE}, {"CZ", CZ}, {"NH1", NH1}, {"NH2", NH2}};
        }
        namespace histidine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 0;
            constexpr int ND1 = 1;
            constexpr int CD2 = 1;
            constexpr int CE1 = 1;
            constexpr int NE2 = 0;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"ND1", ND1}, 
                {"CD2", CD2}, {"CE1", CE1}, {"NE2", NE2}};
        }
        namespace methionine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int SD = 0;
            constexpr int CE = 3;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"SD", SD}, 
                {"CE", CE}};
        }
        namespace cysteine {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int SG = 1;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"SG", SG}};
        }
        namespace proline {
            constexpr int N = 1;
            constexpr int CA = 1;
            constexpr int C = 0;
            constexpr int O = 0;
            constexpr int OXT = 1;
            constexpr int CB = 2;
            constexpr int CG = 2;
            constexpr int CD = 2;
            const std::map<std::string, int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}};
        }
        namespace myristic_acid {
            constexpr int C1 = 0;
            constexpr int O1 = 0;
            constexpr int O2 = 1;
            constexpr int C2 = 2;
            constexpr int C3 = 2;
            constexpr int C4 = 2;
            constexpr int C5 = 2;
            constexpr int C6 = 2;
            constexpr int C7 = 2;
            constexpr int C8 = 2;
            constexpr int C9 = 2;
            constexpr int C10 = 2;
            constexpr int C11 = 2;
            constexpr int C12 = 2;
            constexpr int C13 = 2;
            constexpr int C14 = 3;
            const std::map<std::string, int> get = {{"C1", C1}, {"O1", O1}, {"O2", O2}, {"C2", C2}, {"C3", C3}, {"C4", C4}, {"C5", C5}, {"C6", C6}, {"C7", C7},
                {"C8", C8}, {"C9", C9}, {"C10", C10}, {"C11", C11}, {"C12", C12}, {"C13", C13}, {"C14", C14}};
        }

        // get the number of hydrogen atoms attached to an atom of a specific acid. Example: get.at("GLY").at("CA") = 2
        const std::map<std::string, std::map<std::string, int>> get = {{"GLY", glycine::get}, {"ALA", alanine::get}, {"VAL", valine::get}, 
            {"LEU", leucine::get}, {"ILE", isoleucine::get}, {"PHE", phenylalanine::get}, {"TYR", tyrosine::get}, {"TRP", tryptophan::get}, 
            {"ASP", aspartic_acid::get}, {"GLU", glutamic_acid::get}, {"SER", serine::get}, {"THR", threonine::get}, {"ASN", asparagine::get}, 
            {"GLN", glutamine::get}, {"LYS", lysine::get}, {"ARG", arginine::get}, {"HIS", histidine::get}, {"MET", methionine::get}, 
            {"CYS", cysteine::get}, {"PRO", proline::get}, {"HOH", water::get}, {"MYR", myristic_acid::get}, {"", none::get}};
    }
}