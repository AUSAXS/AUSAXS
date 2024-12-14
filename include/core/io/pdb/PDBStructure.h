#pragma once

#include <io/pdb/PDBFwd.h>
#include <io/IOFwd.h>

#include <io/pdb/Terminate.h>
#include <io/pdb/Footer.h>
#include <io/pdb/Header.h>
#include <utility/TypeTraits.h>

#include <string>
#include <vector>
#include <memory>

namespace ausaxs::io::pdb {
    /**
     * @brief An abstract representation of an input PDBStructure. Handles both reading and writing of PDBStructures, and also stores the relevant PDBStructure data. 
     */
    class PDBStructure {
        public:
            PDBStructure();

            /**
             * @brief Construct a new PDBStructure based on two vectors of atoms. 
             * 
             * @param protein_atoms A vector containing the constituent atoms of the molecule. 
             * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
             */
            PDBStructure(const std::vector<PDBAtom>& protein_atoms, const std::vector<PDBWater>& hydration_atoms);

            /**
             * @brief Construct a new PDBStructure based on two vectors of atoms. 
             * 
             * @param protein_atoms A vector containing the constituent atoms of the molecule. 
             * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
             * @param header The header of this PDBStructure. 
             * @param footer The footer of this PDBStructure. 
             * @param terminate The terminate of this PDBStructure. 
             */
            PDBStructure(const std::vector<PDBAtom>& protein_atoms, const std::vector<PDBWater>& hydration_atoms, const Header& header, const Footer& footer, const Terminate& terminate);

            /**
             * @brief Construct a new PDBStructure based on a input molecular data PDBStructure. 
             * 
             * @param file Path to the input PDBStructure. 
             */
            PDBStructure(const io::File& file);

            /**
             * @brief Update the contents of this PDBStructure.
             * @param patoms The new constituent atoms of the molecule. 
             * @param hatoms The new hydration layer. 
             */
            void update(std::vector<PDBAtom>& patoms, std::vector<PDBWater>& hatoms);

            /**
             * @brief Fill this object with data from a given input data PDBStructure. 
             * @param path Path to the input data PDBStructure. 
             */
            void read(const io::File& path);

            /**
             * @brief Write this PDBStructure to disk. 
             * @param path Path to where this object will be written. 
             */
            void write(const io::File& path);

            /**
             * @brief Get the protein atoms contained in this PDBStructure. 
             */
            const std::vector<PDBAtom>& get_atoms() const;

            /**
             * @brief Get the hydration atoms contained in this PDBStructure. 
             */
            const std::vector<PDBWater> get_waters() const;

            /** 
             * @brief Add an Atom record to this PDBStructure. 
             * @param a The Atom record to be added.
             */
            void add(const PDBAtom& a);

            /** 
             * @brief Add an Atom record to this PDBStructure. 
             * @param a The Atom record to be added.
             */
            void add(PDBAtom&& a);

            /** 
             * @brief Add a Hetatom record to this PDBStructure. 
             * @param w The Hetatom record to be added. 
             */
            void add(const PDBWater& w);

            /** 
             * @brief Add a Hetatom record to this PDBStructure. 
             * @param w The Hetatom record to be added. 
             */
            void add(PDBWater&& w);

            /**
             * @brief Add a Terminate record to this PDBStructure. 
             * @param t The Terminate record to be added. 
             */
            void add(const Terminate& t);

            /**
             * @brief Add a header or footer record to this PDBStructure. 
             * 
             * @param type HEADER or FOOTER.
             * @param s The text string to be added. 
             */
            void add(const RecordType& type, const std::string& s);

            /**
             * @brief Internally updates the consistency of the currently stored data. This ensures this object is in a valid
             *        state for printing. 
             */
            void refresh();

            bool operator==(const PDBStructure& rhs) const;

            bool equals_content(const PDBStructure& rhs) const;

            Header header;
            Footer footer;
            Terminate terminate;
            std::vector<PDBAtom> atoms;
            std::vector<PDBWater> waters;

        private:
            /**
             * @brief Construct a Reader appropriate for the PDBStructure format deduced from the input data PDBStructure. 
             * @param path Path to the input data PDBStructure. 
             */
            std::unique_ptr<io::detail::Reader> construct_reader(const io::File& path);
            /**
             * @brief Construct a Writer appropriate for the PDBStructure format deduced from the output save path.
             * @param path Path to where this object will be written. 
             */
            std::unique_ptr<io::detail::Writer> construct_writer(const io::File& path);
    };
    static_assert(supports_nothrow_move_v<PDBStructure>, "PDBStructure must support nothrow move construction and assignment.");
}