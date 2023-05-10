#pragma once

#include <string>
#include <vector>

#include <data/Terminate.h>
#include <data/Header.h>
#include <data/Footer.h>

class Atom;
class Water;
class Reader;
class Writer;
enum class RecordType;
namespace io {
    class File;
    class ExistingFile;
}

/**
 * @brief An abstract representation of an input ProteinFile. Handles both reading and writing of ProteinFiles, and also stores the relevant ProteinFile data. 
 */
class ProteinFile {
    public:
        /**
         * @brief Default constructor.
         */
        ProteinFile();

        /**
         * @brief Copy constructor.
         */
        ProteinFile(const ProteinFile& ProteinFile);

        /**
         * @brief Move constructor.
         */
        ProteinFile(ProteinFile&& ProteinFile) noexcept;

        /**
         * @brief Constructor. 
         *        Construct a new ProteinFile based on two vectors of atoms. 
         * @param protein_atoms A vector containing the constituent atoms of the molecule. 
         * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
         */
        ProteinFile(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms);

        /**
         * @brief Constructor. 
         *        Construct a new ProteinFile based on two vectors of atoms. 
         * @param protein_atoms A vector containing the constituent atoms of the molecule. 
         * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
         * @param header The header of this ProteinFile. 
         * @param footer The footer of this ProteinFile. 
         * @param terminate The terminate of this ProteinFile. 
         */
        ProteinFile(const std::vector<Atom>& protein_atoms, const std::vector<Water>& hydration_atoms, const Header& header, const Footer& footer, const Terminate& terminate);

        /**
         * @brief Constructor.
         *        Construct a new ProteinFile based on a input molecular data ProteinFile. 
         * @param ProteinFilename Path to the input ProteinFile. 
         */
        ProteinFile(const io::ExistingFile& ProteinFilename);

        /**
         * @brief Destructor.
         */
        ~ProteinFile();

        /**
         * @brief Update the contents of this ProteinFile.
         * @param patoms The new constituent atoms of the molecule. 
         * @param hatoms The new hydration layer. 
         */
        void update(std::vector<Atom>& patoms, std::vector<Water>& hatoms);

        /**
         * @brief Fill this object with data from a given input data ProteinFile. 
         * @param path Path to the input data ProteinFile. 
         */
        void read(const io::ExistingFile& path);

        /**
         * @brief Write this ProteinFile to disk. 
         * @param path Path to where this object will be written. 
         */
        void write(const io::File& path);

        /**
         * @brief Get the protein atoms contained in this ProteinFile. 
         */
        const std::vector<Atom>& get_protein_atoms() const;

        /**
         * @brief Get the hydration atoms contained in this ProteinFile. 
         */
        const std::vector<Water> get_hydration_atoms() const;

        /** 
         * @brief Add an Atom record to this ProteinFile. 
         * @param a The Atom record to be added.
         */
        void add(const Atom& a);

        /** 
         * @brief Add a Hetatom record to this ProteinFile. 
         * @param w The Hetatom record to be added. 
         */
        void add(const Water& w);

        /**
         * @brief Add a Terminate record to this ProteinFile. 
         * @param t The Terminate record to be added. 
         */
        void add(const Terminate& t);

        /**
         * @brief Add a header or footer record to this ProteinFile. 
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

        ProteinFile copy() const;

        ProteinFile& operator=(const ProteinFile& rhs);

        ProteinFile& operator=(ProteinFile&& rhs);

        Header header;
        Footer footer;
        Terminate terminate;
        std::vector<Atom> protein_atoms;
        std::vector<Water> hydration_atoms;

    private:
        std::unique_ptr<Reader> reader;
        std::unique_ptr<Writer> writer;

        /**
         * @brief Construct a Reader appropriate for the ProteinFile format deduced from the input data ProteinFile. 
         * @param path Path to the input data ProteinFile. 
         */
        std::unique_ptr<Reader> construct_reader(const io::ExistingFile& path);
        /**
         * @brief Construct a Writer appropriate for the ProteinFile format deduced from the output save path.
         * @param path Path to where this object will be written. 
         */
        std::unique_ptr<Writer> construct_writer(const io::File& path);
};