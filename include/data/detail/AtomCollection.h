#pragma once

#include <data/record/DataRecordFwd.h>
#include <io/IOFwd.h>

#include <data/record/Terminate.h>
#include <data/record/Footer.h>
#include <data/record/Header.h>

#include <string>
#include <vector>
#include <memory>

namespace data::detail {
    /**
     * @brief An abstract representation of an input AtomCollection. Handles both reading and writing of AtomCollections, and also stores the relevant AtomCollection data. 
     */
    class AtomCollection {
        public:
            /**
             * @brief Default constructor.
             */
            AtomCollection();

            /**
             * @brief Copy constructor.
             */
            AtomCollection(const AtomCollection& AtomCollection);

            /**
             * @brief Move constructor.
             */
            AtomCollection(AtomCollection&& AtomCollection) noexcept;

            /**
             * @brief Construct a new AtomCollection based on two vectors of atoms. 
             * 
             * @param protein_atoms A vector containing the constituent atoms of the molecule. 
             * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
             */
            AtomCollection(const std::vector<record::Atom>& protein_atoms, const std::vector<record::Water>& hydration_atoms);

            /**
             * @brief Construct a new AtomCollection based on two vectors of atoms. 
             * 
             * @param protein_atoms A vector containing the constituent atoms of the molecule. 
             * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
             * @param header The header of this AtomCollection. 
             * @param footer The footer of this AtomCollection. 
             * @param terminate The terminate of this AtomCollection. 
             */
            AtomCollection(const std::vector<record::Atom>& protein_atoms, const std::vector<record::Water>& hydration_atoms, const record::Header& header, const record::Footer& footer, const record::Terminate& terminate);

            /**
             * @brief Construct a new AtomCollection based on a input molecular data AtomCollection. 
             * 
             * @param file Path to the input AtomCollection. 
             */
            AtomCollection(const io::File& file);

            /**
             * @brief Destructor.
             */
            ~AtomCollection();

            /**
             * @brief Update the contents of this AtomCollection.
             * @param patoms The new constituent atoms of the molecule. 
             * @param hatoms The new hydration layer. 
             */
            void update(std::vector<record::Atom>& patoms, std::vector<record::Water>& hatoms);

            /**
             * @brief Fill this object with data from a given input data AtomCollection. 
             * @param path Path to the input data AtomCollection. 
             */
            void read(const io::File& path);

            /**
             * @brief Write this AtomCollection to disk. 
             * @param path Path to where this object will be written. 
             */
            void write(const io::File& path);

            /**
             * @brief Get the protein atoms contained in this AtomCollection. 
             */
            const std::vector<record::Atom>& get_protein_atoms() const;

            /**
             * @brief Get the hydration atoms contained in this AtomCollection. 
             */
            const std::vector<record::Water> get_hydration_atoms() const;

            /** 
             * @brief Add an Atom record to this AtomCollection. 
             * @param a The Atom record to be added.
             */
            void add(const record::Atom& a);

            /** 
             * @brief Add an Atom record to this AtomCollection. 
             * @param a The Atom record to be added.
             */
            void add(record::Atom&& a);

            /** 
             * @brief Add a Hetatom record to this AtomCollection. 
             * @param w The Hetatom record to be added. 
             */
            void add(const record::Water& w);

            /** 
             * @brief Add a Hetatom record to this AtomCollection. 
             * @param w The Hetatom record to be added. 
             */
            void add(record::Water&& w);

            /**
             * @brief Add a Terminate record to this AtomCollection. 
             * @param t The Terminate record to be added. 
             */
            void add(const record::Terminate& t);

            /**
             * @brief Add a header or footer record to this AtomCollection. 
             * 
             * @param type HEADER or FOOTER.
             * @param s The text string to be added. 
             */
            void add(const record::RecordType& type, const std::string& s);

            /**
             * @brief Internally updates the consistency of the currently stored data. This ensures this object is in a valid
             *        state for printing. 
             */
            void refresh();

            AtomCollection copy() const;

            AtomCollection& operator=(const AtomCollection& rhs);

            AtomCollection& operator=(AtomCollection&& rhs);

            bool operator==(const AtomCollection& rhs) const;

            bool equals_content(const AtomCollection& rhs) const;

            record::Header header;
            record::Footer footer;
            record::Terminate terminate;
            std::vector<record::Atom> protein_atoms;
            std::vector<record::Water> hydration_atoms;

        private:
            std::unique_ptr<io::detail::Reader> reader;
            std::unique_ptr<io::detail::Writer> writer;

            /**
             * @brief Construct a Reader appropriate for the AtomCollection format deduced from the input data AtomCollection. 
             * @param path Path to the input data AtomCollection. 
             */
            std::unique_ptr<io::detail::Reader> construct_reader(const io::File& path);
            /**
             * @brief Construct a Writer appropriate for the AtomCollection format deduced from the output save path.
             * @param path Path to where this object will be written. 
             */
            std::unique_ptr<io::detail::Writer> construct_writer(const io::File& path);
    };
}