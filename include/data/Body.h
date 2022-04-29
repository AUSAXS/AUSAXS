#pragma once

#include <vector>
#include <map>
#include <memory>

#include <data/Atom.h>
#include <hydrate/Grid.h>
#include <io/File.h>
#include <utility/Constants.h>
#include <histogram/ScatteringHistogram.h>
#include <data/StateManager.h>

class Body {
  public:
    /**
     * @brief Default constructor.
     * 
     * This is not a very efficient way of constructing a Body.
     */
    Body();

    /** 
     * @brief Create a new collection of atoms (body) from the input .pdb or .xml file. 
     * 
     * @param path path to the input file. 
     * @param signaller a signalling object to signal changes of state
     */
    explicit Body(std::string path);

    /**
     * @brief Create a new collection of atoms (body) based on two vectors
     */
    explicit Body(const std::vector<Atom>& protein_atoms, const std::vector<Hetatom>& hydration_atoms = {});

    /**
     * @brief Copy constructor. 
     */
    Body(const Body& body);

    /**
     * @brief Move constructor. 
     */
    Body(Body&& body);

    ~Body();

    /** 
     * @brief Writes this body to disk.
     * 
     * @param path path to the destination. 
     */
    void save(std::string path);

    /**
     * @brief Get the distances between each atom.
     */
    std::shared_ptr<ScatteringHistogram> get_histogram();

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this body. Note that the previous one will be deleted.
     */
    void generate_new_hydration();

    /**
     * @brief Get a reference to the constituent atoms.
     */
    const std::vector<Atom>& get_protein_atoms() {return protein_atoms;}

    /**
     * @brief Get a referece to the hydration atoms.
     */
    const std::vector<Hetatom>& get_hydration_atoms() {return hydration_atoms;}

    /** 
     * @brief Calculate the center-mass coordinates for the body.
     * @return The center-mass (x, y, z) coordinates. 
     */
    Vector3 get_cm() const;

    /**
     * @brief Calculate the volume of this body based on its constituent amino acids
     */
    double get_volume_acids() const;

    /**
     * @brief Calculate the volume of this body based on the number of grid bins it spans
     */
    double get_volume_grid();

    /**
     * @brief Calculate the volume of this body based on the number of C-alpha atoms
     */
    double get_volume_calpha() const;

    /**
     * @brief Get the grid representation of this body. 
     */
    std::shared_ptr<Grid> get_grid();

    /**
     * @brief Generate a PDB file at @p path showing the filled grid volume.
     */
    void generate_volume_file(std::string path);

    /**
     * @brief Calculate the total mass of this body in Daltons.
     */
    double get_mass() const;

    /**
     * @brief Create a grid and fill it with the atoms of this body. 
     */
    std::shared_ptr<Grid> create_grid();

    /**
     * @brief Center this Body on origo. 
     */
    void center();

    /** 
     * @brief Move the entire body by a vector.
     * @param v the translation vector
     */
    void translate(const Vector3& v);

    /**
     * @brief Rotate all atoms by a given rotation matrix.
     * 
     * @param R The rotation matrix. 
     */
    void rotate(const Matrix<double>& R);
    
    /**
     * @brief Rotate all atoms @a rad radians about the axis @a axis. 
     * 
     * @param axis the rotation axis. 
     * @param rad the amount to rotate in radians. 
     */
    void rotate(const Vector3& axis, double rad);

    /**
     * ! Not implemented
     * @brief Euler angle rotation of all atoms. 
     * 
     * @param alpha radians to rotate about the z-axis.
     * @param beta radians to rotate about the y-axis. 
     * @param gamma radians to rotate about the x-axis. 
     */
    void rotate(double alpha, double beta, double gamma);

    /** 
     * @brief Calculate the distances between each pair of atoms. 
     */
    void calc_histogram();

    /**
     * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
     */
    void update_effective_charge();

    /**
     * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
     * 
     * @param charge the charge to be subtracted.
     */
    void update_effective_charge(double charge);

    /**
     * @brief Register a probe (listener) to this object, which will be notified of state changes. 
     */
    void register_probe(std::shared_ptr<StateManager::BoundSignaller> signal);

    /**
     * @brief Assign another body to this object. 
     */
    Body& operator=(const Body& rhs);

    /**
     * @brief Check if this object is equal to another. 
     */
    bool operator==(const Body& rhs) const;

    /**
     * @brief Get the File backing this object. 
     */
    std::shared_ptr<File> get_file() const {return file;}

    /**
     * @brief Signal that this object has changed its external state.
     *        This triggers recalculating all external distances between this body and everything else the next time a histogram is requested. 
     */
    void changed_external_state() const;

    /**
     * @brief Signal that this object has changed its internal state.
     *        This triggers recalculating all distances, both external and internal, between this body and everything else the next time a histogram is requested. 
     */
    void changed_internal_state() const;    

    std::shared_ptr<StateManager::Signaller> signal = std::make_shared<StateManager::UnboundSignaller>(); 
  private:
    std::shared_ptr<File> file = nullptr;                     // The file backing this body
    std::shared_ptr<Grid> grid = nullptr;                     // The grid representation of this body
    std::shared_ptr<ScatteringHistogram> histogram = nullptr; // An object representing the distances between atoms

    // The signalling object to signal a change of state. The default doesn't do anything, and must be overriden by a proper Signaller object.  
    // shared_ptr<StateManager::Signaller> signal = std::make_shared<StateManager::UnboundSignaller>(); 

  public: 
    size_t uid;                           // An unique identifier for this body
    std::vector<Atom>& protein_atoms;          // Atoms of the body itself
    std::vector<Hetatom>& hydration_atoms;     // Hydration layer
    bool updated_charge = false;          // True if the effective charge of each atom has been updated to reflect the volume they occupy, false otherwise
    bool centered = false;                // True if this object is centered, false otherwise
    inline static size_t uid_counter = 0; // The unique counter. 
};