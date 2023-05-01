#include <utility/Constants.h>
#include <utility/StringUtils.h>

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

    namespace volume {
        const saxs::detail::SimpleMap<double> amino_acids = std::unordered_map<std::string, double>{
            {"GLY", 66.4}, {"ALA", 91.5}, {"VAL", 141.7}, {"LEU", 167.9}, {"ILE", 168.8}, {"PHE", 203.5}, {"TYR", 203.6}, {"TRP", 237.6}, 
            {"ASP", 113.6}, {"GLU", 140.6}, {"SER", 99.1}, {"THR", 122.1}, {"ASN", 135.2}, {"GLN", 161.1}, {"LYS", 176.2}, {"ARG", 180.8}, 
            {"HIS", 167.3}, {"MET", 170.8}, {"CYS", 105.6}, {"PRO", 129.3}
        };
    }

    namespace mass {
        // The fake element "M" is for compatibility with the GROMACS tip4p water model. 
        const saxs::detail::SimpleMap<double> atomic = std::unordered_map<std::string, double>{
            {"H", 1.0079}, {"He", 4.0026}, {"Li", 6.941}, {"Be", 9.0122}, {"B", 10.811}, {"C", 12.0107}, {"N", 14.0067}, {"O", 15.9994}, {"F", 18.9984}, 
            {"Ne", 20.1797}, {"Na", 22.9897}, {"Mg", 24.305}, {"Al", 26.9815}, {"Si", 28.0855}, {"P", 30.9738}, {"S", 32.065}, {"Cl", 35.453}, 
            {"Ar", 39.948}, {"K", 39.0983}, {"Ca", 40.078}, {"Sc", 44.9559}, {"Ti", 47.867}, {"V", 50.9415}, {"Cr", 51.9961}, {"Mn", 54.938}, 
            {"Fe", 55.845}, {"Co", 58.9332}, {"Ni", 58.6934}, {"Cu", 63.546}, {"Zn", 65.39}, {"W", 183.84},
            {"M", 0}
        };
    }

    namespace charge {
        const saxs::detail::SimpleMap<unsigned int> atomic = std::unordered_map<std::string, unsigned int>{
            {"e", 1}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, 
            {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, 
            {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"W", 74}, 
            {"M", 0}
        };
    }

    namespace valence {
        const saxs::detail::SimpleMap<unsigned int> atomic = std::unordered_map<std::string, unsigned int>{
            {"H", 1}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"Ne", 0}, {"S", 2}, {"P", 1}, {"Cl", 1}, {"M", 0}
        };
    }

    namespace hydrogen_atoms {
        parser::residue::ResidueStorage residues;
    }

    namespace symbols {
        std::string hydrogen = "H";
        std::string carbon = "C";
        std::string nitrogen = "N";
        std::string oxygen = "O";
    }
}