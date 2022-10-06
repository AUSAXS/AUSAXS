#pragma once

#include <vector>
#include <string>
#include <map>

#include <utility/ResidueMap.h>

namespace parser {
    namespace residue {
        namespace detail {
            struct Atom {
                Atom(std::string name, std::string altname, std::string symbol);

                Atom(std::string name, int charge);

                void add_bond(std::string symbol, unsigned int order);

                std::string to_string() const;

                std::string name, altname, symbol;
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

                    void add_atom(std::string name, std::string altname, std::string symbol);

                    void add_atom(std::string name, int charge);

                    void apply_bond(const std::vector<Bond>& bonds);

                    void apply_bond(const Bond& bond);

                    std::string to_string() const;

                    saxs::detail::ResidueMap to_map() const;

                    static Residue parse(std::string filename);

                private: 
                    std::string name;
                    std::map<std::string, unsigned int> name_map;
                    std::vector<Atom> atoms;        
            };

            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Atom& a);
            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Bond& b);
            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Residue& l);
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
                saxs::detail::ResidueMap& get(std::string name);

            private: 
                /**
                 * @brief Insert a residue into the storage. 
                 */
                void insert(std::string name, saxs::detail::ResidueMap residue);

                /**
                 * @brief Initialize this storage. All residue files present in the storage directory will be loaded. 
                 */
                void initialize();

                /**
                 * @brief Download and load a residue from the web. 
                 */
                void download_residue(std::string name);

                /**
                 * @brief Write a residue to the auto-loaded file. 
                 */
                void write_residue(std::string name);

                std::map<std::string, saxs::detail::ResidueMap> data;
        };
    }
}