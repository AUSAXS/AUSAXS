#pragma once

#include <utility/ResidueMap.h>

#include <vector>
#include <string>
#include <unordered_map>

namespace constants {enum class atom_t;}
namespace io {class ExistingFile;}
namespace parser {
    namespace residue {
        namespace detail {
            /**
             * @brief A simple representation of an atom.
             */
            struct Atom {
                /**
                 * @brief Constructor.
                 * 
                 * @param name The name of the atom.
                 * @param altname The alternate name of the atom.
                 * @param atom The atom type. This is necessary to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
                 */
                Atom(const std::string& name, const std::string& altname, constants::atom_t symbol);

                /**
                 * @brief Constructor.
                 * 
                 * @param name The name of the atom.
                 * @param charge The charge of the atom.
                 * @param atom The atom type. This is necessary to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
                 */
                Atom(const std::string& name, int charge, constants::atom_t symbol);

                /**
                 * @brief Add a bond from this atom to another.
                 * 
                 * @param symbol The type of atom this instance is bonded to.
                 * @param order The order of the bond.
                 */
                void add_bond(constants::atom_t symbol, unsigned int order);
                
                std::string to_string() const;

                std::string name, altname;
                constants::atom_t atom;
                unsigned int valency;
                unsigned int hydrogen_bonds = 0;
            };

            /**
             * @brief A simple representation of a bond. 
             */
            struct Bond {
                /**
                 * @brief Constructor.
                 * 
                 * @param name1 the name of the first atom.
                 * @param name2 the name of the second atom.
                 * @param order the order of the bond.
                 */
                Bond(const std::string& name1, const std::string& name2, unsigned int order);

                std::string to_string() const;

                /**
                 * @brief Parse an order symbol.
                 */
                static unsigned int parse_order(const std::string& order);

                std::string name1, name2;
                unsigned int order;
            };

            /**
             * @brief A simple representation of a residue.
             */
            class Residue {
                public: 
                    /**
                     * @brief Construct a new residue with the given name.
                     */
                    Residue(const std::string& name);

                    /**
                     * @brief Constructor.
                     * 
                     * @param name The name of the residue.
                     * @param atoms The atoms in the residue.
                     * @param bonds The bonds in the residue.
                     */
                    Residue(const std::string& name, std::vector<Atom> atoms, std::vector<Bond> bonds);

                    /**
                     * @brief Add an atom to the residue.
                     * 
                     * @param name The name of the atom.
                     * @param altname The alternate name of the atom.
                     * @param atom The atom type.
                     */
                    void add_atom(const std::string& name, const std::string& altname, constants::atom_t atom);

                    /**
                     * @brief Add an atom to the residue.
                     * 
                     * @param name The name of the atom.
                     * @param charge The charge of the atom.
                     * @param atom The atom type.
                     */
                    void add_atom(const std::string& name, int charge, constants::atom_t atom);

                    /**
                     * @brief Add a bond to the residue.
                     * 
                     * @param name1 the name of the first atom.
                     * @param name2 the name of the second atom.
                     * @param order the order of the bond.
                     */
                    void apply_bond(const std::vector<Bond>& bonds);

                    /**
                     * @brief Add a bond to the residue.
                     * 
                     * @param name1 the name of the first atom.
                     * @param name2 the name of the second atom.
                     * @param order the order of the bond.
                     */
                    void apply_bond(const Bond& bond);

                    /**
                     * @brief Parse a residue from a file.
                     */
                    saxs::detail::ResidueMap to_map() const;

                    /**
                     * @brief Parse a residue from a file.
                     */
                    static Residue parse(const io::ExistingFile& filename);

                    std::string to_string() const;

                private: 
                    std::string name;
                    std::unordered_map<std::string, int> name_map;
                    std::vector<Atom> atoms;        
            };

            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Atom& a);
            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Bond& b);
            std::ostream& operator<<(std::ostream& os, const parser::residue::detail::Residue& l);
        }

        /**
         * @brief A storage container for residues.
         */
        class ResidueStorage {
            public: 
                /**
                 * @brief Default constructor. The storage will be initialized with all residue files present in the storage directory.
                 */
                ResidueStorage();

                ~ResidueStorage();

                /**
                 * @brief Get a residue from the storage. If the residue is not found, it will be downloaded. 
                 */
                saxs::detail::ResidueMap& get(const std::string& name);

            private: 
                /**
                 * @brief Insert a residue into the storage. 
                 */
                void insert(const std::string& name, const saxs::detail::ResidueMap& residue);

                /**
                 * @brief Initialize this storage. All residue files present in the storage directory will be loaded. 
                 */
                void initialize();

                /**
                 * @brief Download and load a residue from the web. 
                 */
                void download_residue(const std::string& name);

                /**
                 * @brief Write a residue to the auto-loaded file. 
                 */
                void write_residue(const std::string& name);

                std::unordered_map<std::string, saxs::detail::ResidueMap> data;
        };
    }
}