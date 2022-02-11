#pragma once

#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include "data/Record.h"
#include "data/Terminate.h"
#include "data/Header.h"
#include "data/Footer.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "io/Reader.h"
#include "io/Writer.h"
#include "io/PDBWriter.h"
#include "io/PDBReader.h"

using std::vector, std::string, std::unique_ptr, std::shared_ptr; 

/**
 * @brief \class File
 * An abstract representation of an input file. Handles both reading and writing of files, and also stores the relevant file data. 
 */
class File {
  public: 
    /**
     * @brief Default constructor.
     */
    File() {}

    /**
     * @brief Copy constructor.
     */
    File(const File& file);

    /**
     * @brief Move constructor.
     */
    File(const File&& file) noexcept;

    /**
     * @brief Constructor. 
     *        Construct a new File based on two vectors of atoms. 
     * @param protein_atoms A vector containing the constituent atoms of the molecule. 
     * @param hydration_atoms A vector containing the water molecules for an existing hydration layer. 
     */
    File(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms);

    /**
     * @brief Constructor.
     *        Construct a new File based on a input molecular data file. 
     * @param filename Path to the input file. 
     */
    File(string filename);

    /**
     * @brief Destructor.
     */
    ~File();

    /**
     * @brief Update the contents of this file.
     * @param patoms The new constituent atoms of the molecule. 
     * @param hatoms The new hydration layer. 
     */
    void update(vector<Atom>& patoms, vector<Hetatom>& hatoms);

    /**
     * @brief Fill this object with data from a given input data file. 
     * @param path Path to the input data file. 
     */
    void read(const string& path);

    /**
     * @brief Write this file to disk. 
     * @param path Path to where this object will be written. 
     */
    void write(const string path);

    /**
     * @brief Get the protein atoms contained in this File. 
     */
    const vector<Atom>& get_protein_atoms() const;

    /**
     * @brief Get the hydration atoms contained in this File. 
     */
    const vector<Hetatom> get_hydration_atoms() const;

    /** 
     * @brief Add an Atom record to this file. 
     * @param r The Atom record to be added.
     */
    void add(const Atom r);

    /** 
     * @brief Add a Hetatom record to this file. 
     * @param r The Hetatom record to be added. 
     */
    void add(const Hetatom r);

    /**
     * @brief Add a Terminate record to this file. 
     * @param r The Terminate record to be added. 
     */
    void add(const Terminate);

    /**
     * @brief Add a header or footer record to this file. 
     * @param type HEADER or FOOTER.
     * @param s The text string to be added. 
     */
    void add(const string type, const string s);

    /**
     * @brief Internally updates the consistency of the currently stored data. This ensures this object is in a valid
     *        state for printing. 
     */
    void refresh();

    File copy() const;

    File& operator=(const File& rhs);

    Header header;
    Footer footer;
    Terminate terminate;
    vector<Atom> protein_atoms;
    vector<Hetatom> hydration_atoms;

  private:
    std::unique_ptr<Reader> reader;
    std::unique_ptr<Writer> writer;

    /**
     * @brief Construct a Reader appropriate for the file format deduced from the input data file. 
     * @param path Path to the input data file. 
     */
    std::unique_ptr<Reader> construct_reader(const string& path);
    /**
     * @brief Construct a Writer appropriate for the file format deduced from the output save path.
     * @param path Path to where this object will be written. 
     */
    std::unique_ptr<Writer> construct_writer(const string& path);
};