#pragma once

#include "data/Atom.h"
#include "PlacementStrategy.h"
#include "CullingStrategy.h"
#include "settings.h"
#include "Exceptions.h"
#include "hydrate/GridMember.h"

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
    vector<GridMember<T>> add(const vector<T>& atoms) {
        static_assert(std::is_base_of<Atom, T>::value, "Argument type must be derivable from Atom!");
        vector<GridMember<T>> v(atoms.size());
        size_t index = 0;
        for (auto const& a : atoms) {
            v[index++] = add(a);
        }
        return v;
    }

    /** 
     * @brief Add a single atom to the grid. 
     *        This is a constant-time operation. 
     */
    GridMember<Atom> add(const Atom& atom, const bool& expand = false);

    /** 
     * @brief Add a single hetatom to the grid. 
     *        This is a constant-time operation. 
     */
    GridMember<Hetatom> add(const Hetatom& atom, const bool& expand = false);

    /**
     * @brief Remove a single atom from the grid.
     *        This is linear in the number of stored elements. 
     */
    void remove(const Atom& atom);

    /**
     * @brief Remove a single water molecule from the grid.
     *        This is linear in the number of stored elements. 
     */
    void remove(const Hetatom& atom);

    /**
     * @brief Remove multiple water molecule from the grid.
     *        This is linear in the number of stored elements times the size of the input vector. 
     */
    void remove(const vector<Hetatom>& atom);

    /** 
     * @brief Expand all member atoms and water molecules into actual spheres based on the radii ra and rh. 
     */
    void expand_volume();

    /** 
     * @brief Expand a single member atom into an actual sphere.
     */
    void expand_volume(GridMember<Atom>& atom);

    /** 
     * @brief Expand a single member atom into an actual sphere.
     */
    void expand_volume(GridMember<Hetatom>& atom);

    /** 
     * @brief Deflate all member atoms and water molecules into actual spheres based on the radii ra and rh. 
     */
    void deflate_volume();

    /** 
     * @brief Deflate a single member atom into an actual sphere.
     */
    void deflate_volume(GridMember<Atom>& atom);

    /** 
     * @brief Deflate a single member atom into an actual sphere.
     */
    void deflate_volume(GridMember<Hetatom>& atom);

    /**
     * @brief Generate a new hydration layer for the grid.
     * @return Pointers to the new water molecules. 
     */
    vector<Hetatom> hydrate();

    /**
     * @brief Identify possible hydration binding locations for the structure. 
     * @return A list of possible (binx, biny, binz) locations.
     */
    vector<GridMember<Hetatom>> find_free_locs();

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
     * @brief Get all water molecules from this grid. 
     */
    vector<Hetatom> get_waters() const;

    /**
     * @brief Get all atoms from this grid. 
     */
    vector<Atom> get_atoms() const;

    /**
     * @brief Get the total volume spanned by the atoms in this grid in Å^3.
     *        Water molecules are ignored.  
     */
    double get_volume();

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

    // We define our own comparison method for the map. We must do this since it must be able to hold duplicate atoms (water molecules),
    // so we instead compare their unique identifiers
    class Comparator {public: bool operator() (const Atom& l, const Atom& r) const {return l.uid < r.uid;}};

    struct MapVal {
        vector<int> loc; // the bin location of the Atom key
        bool expanded_volume; // whether the volume of this location has been expanded

        // the two operator overloads makes this struct act just like a vector, but with an additional bool available when needed. 
        int operator[](int index) {return loc[index];}
        int operator[](int index) const {return loc[index];}
    };

    vector<vector<vector<char>>> grid; // The actual grid. Datatype is char since we need at least four different values.
    std::list<GridMember<Atom>> a_members; // A list of all member atoms and where they are located.
    std::list<GridMember<Hetatom>> w_members; // A list of all member water molecules and where they are located. 
    int volume = 0; // The number of bins covered by the members, i.e. the actual volume in the unit (width)^3
    int ra = 0; // Radius of each atom represented as a number of bins
    int rh = 0; // Radius of each water molecule represented as a number of bins

private:
    Vector3 base; // base point of this grid
    double width; // distance between each grid point
    vector<int> bins; // the number of bins in each dimension
    unique_ptr<PlacementStrategy> water_placer; // the strategy for placing water molecules
    unique_ptr<CullingStrategy> water_culler; // the strategy for culling the placed water molecules

    /** 
     * @brief Expand a single member atom into an actual sphere.
     * @param loc The bin location of the atom. 
     * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
     */
    void expand_volume(const vector<int>& val, const bool is_water);

    /** 
     * @brief Deflate a single member atom into an actual sphere.
     * @param loc The bin location of the atom. 
     * @param is_water If the atom is a water molecule. Used to determine which marker to use in the grid. 
     */
    void deflate_volume(const vector<int>& loc, const bool is_water);
};