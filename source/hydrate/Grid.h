#pragma once

// includes
#include "data/Atom.h"
#include "PlacementStrategy.h"
#include "CullingStrategy.h"
#include "settings.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Grid {
public:
    /**
     * @brief Construct a new Grid object with standard atomic radii.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     */
    Grid(Vector3 base, double width, int bins) : Grid(base, width, {bins, bins, bins}, setting::grid::ra, setting::grid::rh, setting::grid::psc, setting::grid::csc) {};

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     * @param radius the radius of each atom.
     */
    Grid(Vector3 base, double width, int bins, int radius) : Grid(base, width, {bins, bins, bins}, radius, radius, setting::grid::psc, setting::grid::csc) {};

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in each dimension. 
     * @param ra the radius of each atom.
     * @param rh the radius of each water molecule.
     */
    Grid(Vector3 base, double width, vector<int> bins, double ra, double rh, setting::grid::PlacementStrategyChoice psc, setting::grid::CullingStrategyChoice csc);

    /** 
     * @brief Add a set of atoms to the grid. 
     * @param atoms the set of atoms to add to this grid.
     */
    template<typename T>
    void add(const vector<shared_ptr<T>>& atoms) {
        static_assert(std::is_base_of<Atom, T>::value, "Argument type must be derivable from Atom!");
        for (auto const& a : atoms) {
            add(a);
        }
    }

    /** 
     * @brief Add a single atom to the grid. 
     * @param atom the atom to be added. 
     */
    void add(const shared_ptr<Atom> atom);

    /**
     * @brief Remove a single atom from the grid.
     * @param atom the atom to be removed.
     */
    void remove(const shared_ptr<Atom> atom);

    /** 
     * @brief Expand the member atoms into actual spheres based on the radii ra and rh. 
     */
    void expand_volume();

    /** 
     * @brief Expand a single member atom into an actual sphere.
     * @param atom the member atom to be expanded. 
     */
    void expand_volume(const shared_ptr<Atom> atom);

    /**
     * @brief Generate a new hydration layer for the grid.
     * @return Pointers to the new water molecules. 
     */
    vector<shared_ptr<Hetatom>> hydrate();

    /**
     * @brief Identify possible hydration binding locations for the structure. 
     * @return A list of possible (binx, biny, binz) locations.
     */
    vector<shared_ptr<Hetatom>> find_free_locs();

    /**
     * @brief Create the smallest possible box containing the center points of all member atoms.
     * @return vector<vector<int>> An index pair (min, max) for each dimension (shape: [3][2]). 
     */
    vector<vector<int>> bounding_box() const;

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
    vector<shared_ptr<Atom>> get_hydration_atoms() const;

    /**
     * @brief Get all protein atoms from this grid. 
     */
    vector<shared_ptr<Atom>> get_protein_atoms() const;

    /**
     * @brief Get the total volume spanned by the atoms in this grid in Å^3.
     *        Water molecules are ignored.  
     */
    double get_volume() {
        if (!vol_expanded) {expand_volume();} 
        return pow(width, 3)*volume;
    }

    /**
     * @brief Get the number of bins in each dimension.
     */
    const vector<int> get_bins() const {return bins;}

    /**
     * @brief Get the width of each bin.
     */
    double get_width() const {return width;}

    /**
     * @brief Convert a vector of bin locations (binx, biny, binz) to a vector of absolute coordinates (x, y, z).
     * @param v the bin location.
     * @return The xyz location.
     */
    Vector3 to_xyz(const vector<int>& v) const;

    /**
     * @brief Convert a vector of absolute coordinates (x, y, z) to a vector of bin locations.
     * @param v the xyz location.
     * @return The bin location. 
     */
    vector<int> to_bins(const Vector3& v) const;

    class Comparator {
        public: bool operator() (const shared_ptr<Atom>& l, const shared_ptr<Atom>& r) const {return *l < *r;}
    };

    vector<vector<vector<char>>> grid; // the actual grid. Datatype is char since we need at least four different values
    std::map<const shared_ptr<Atom>, vector<int>> members; // a map of all members of this grid and where they are located
    int volume = 0; // the number of bins covered by the members, i.e. the actual volume in the unit (width)^3
    int ra = 0; // radius of each atom represented as a number of bins
    int rh = 0; // radius of each water molecule represented as a number of bins

private:
    Vector3 base; // base point of this grid
    double width; // distance between each grid point
    vector<int> bins; // the number of bins in each dimension
    bool vol_expanded = false; // a flag determining if the volume has been expanded 
    unique_ptr<PlacementStrategy> water_placer; // the strategy for placing water molecules
    unique_ptr<CullingStrategy> water_culler; // the strategy for culling the placed water molecules
};