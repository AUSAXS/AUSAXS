#pragma once

// includes
#include <TVector3.h>
#include "data/Atom.cpp"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class Grid {
public:
    /**
     * @brief Construct a new Grid object with standard atomic radii.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     */
    Grid(TVector3 base, double width, int bins) : Grid(base, width, {bins, bins, bins}, sqrt(8), 1) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     * @param radius the radius of each atom.
     */
    Grid(TVector3 base, double width, int bins, int radius) : Grid(base, width, {bins, bins, bins}, radius, radius) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in each dimension. 
     * @param ra the radius of each atom.
     * @param rh the radius of each water molecule.
     */
    Grid(TVector3 base, double width, vector<int> bins, double ra, double rh);

    /** 
     * @brief Add a set of atoms to the grid. 
     * @param atoms the set of atoms to add to this grid.
     */
    void add(vector<shared_ptr<Atom>>* atoms);

    /** 
     * @brief Expand the member atoms into actual spheres based on the radii ra and rh. 
     */
    void expand_volume();

    /**
     * @brief Generate a new hydration layer for the grid.  
     * @param reduce reduces the number of generated HOH molecules by this factor. Use 0 for no reduction. 
     * @return Pointers to the new water molecules. 
     */
    vector<shared_ptr<Atom>> hydrate(int reduce);

    /**
     * @brief Identify possible hydration binding locations for the structure. 
     * @return A list of possible (binx, biny, binz) locations.
     */
    vector<vector<int>> find_free_locs();

    /**
     * @brief Create the smallest possible box containing the center points of all member atoms.
     * @return vector<vector<int>> An index pair (min, max) for each dimension (shape: [3][2]). 
     */
    vector<vector<int>> bounding_box() const;

    /** 
     * @brief Get a pointer to the internal grid. 
     */
    vector<vector<vector<char>>>* get_grid();

    /**
     * @brief Set the radius of all atoms (not water molecules!).
     * @param radius The new radius in Ångström.
     */
    void set_radius_atoms(double radius);

    /**
     * @brief Set the radius for all water molecules.
     * @param radius The new radius in Ångström.
     */
    void set_radius_water(double radius);

    /**
     * @brief Get all hydration atoms from this grid. 
     */
    vector<Atom*> get_hydration_atoms() const;

    /**
     * @brief Get the total volume spanned by the atoms in this grid. 
     */
    double get_volume();

protected:
    TVector3 base; // base point of this grid
    double width; // distance between each grid point
    vector<vector<vector<char>>> grid; // the actual grid. Datatype is char since we need at least four different values
    std::map<Atom, vector<int>> members; // a map of all members of this grid and where they are located
    vector<int> bins; // the number of bins in each dimension
    int volume = 0; // the number of bins covered by the members, i.e. the actual volume in the unit (width)^3
    int ra = 0; // radius of each atom represented as a number of bins
    int rh = 0; // radius of each water molecule represented as a number of bins
    bool vol_expanded = false; // a flag determining if the volume has been expanded 

private:
    /**
     * @brief Check if a water molecule can be placed at the given location. 
     *        This checks collisions with both other water molecules and other atoms. 
     * @param loc the location to be checked. 
     * @param other_molecules the other water molecules which have already been placed.
     * @return True if this is an acceptable location, false otherwise.
     */
    bool check_collisions(const vector<int> loc, const vector<vector<int>> other_molecules) const;

    /** 
     * @brief Expand a single member atom into an actual sphere.
     * @param atom the member atom to be expanded. 
     */
    void expand_volume(const Atom atom);

    /** 
     * @brief Add a single atom to the grid. 
     * @param atoms the atom to be added. 
     */
    void add(shared_ptr<Atom> atom);

    /**
     * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
     * @param v the xyz location.
     * @return The bin location. 
     */
    vector<int> to_bins(TVector3 v);

    /**
     * @brief Convert a vector of bin locations (binx, biny, binz) to a vector of absolute coordinates (x, y, z).
     * @param v the bin location.
     * @return The xyz location.
     */
    TVector3 to_xyz(vector<int> v);
};