#pragma once

#include <vector>
#include <string>
#include <map>

#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <curl/curl.h>

namespace parser {
    namespace ligand {
        namespace detail {
            struct Atom {
                Atom(std::string name, std::string symbol);

                void add_bond(std::string symbol, unsigned int order);

                std::string to_string() const;

                std::string name, symbol;
                unsigned int valency;
                unsigned int hydrogen_bonds = 0;
            };

            struct Bond {
                Bond(std::string name1, std::string name2, unsigned int order);

                std::string to_string() const;

                static unsigned int parse_order(std::string order);

                std::string name1, name2;
                unsigned int order;
            };

            class Residue {
                public: 
                    Residue(std::string name);

                    Residue(std::string name, std::vector<Atom> atoms, std::vector<Bond> bonds);

                    void add_atom(std::string name, std::string symbol);

                    void apply_bond(const std::vector<Bond>& bonds);

                    void apply_bond(const Bond& bond);

                    std::string to_string() const;

                    std::map<std::string, unsigned int> to_map() const;

                    static Residue parse(std::string filename);

                private: 
                    std::string name;
                    std::map<std::string, unsigned int> name_map;
                    std::vector<Atom> atoms;        
            };

            std::ostream& operator<<(std::ostream& os, const parser::ligand::detail::Atom& a);
            std::ostream& operator<<(std::ostream& os, const parser::ligand::detail::Bond& b);
            std::ostream& operator<<(std::ostream& os, const parser::ligand::detail::Residue& l);
        }

        class ResidueStorage {
            public: 
                /**
                 * @brief Default constructor. The storage will be initialized with all residue files present in the storage directory.
                 */
                ResidueStorage();

                /**
                 * @brief Get a residue from the storage. If the residue is not found, it will be downloaded. 
                 */
                std::map<std::string, unsigned int>& get(std::string name);

            private: 
                /**
                 * @brief Insert a residue into the storage. 
                 */
                void insert(std::string name, std::map<std::string, unsigned int> ligand);

                /**
                 * @brief Initialize this storage. All residue files present in the storage directory will be loaded. 
                 */
                void initialize();

                /**
                 * @brief Download and load a residue from the web. 
                 */
                void download_residue(std::string name);

                std::map<std::string, std::map<std::string, unsigned int>> data;
        };
    }
}