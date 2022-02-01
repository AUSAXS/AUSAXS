#pragma once

#include <list>

#include "data/Atom.h"
#include "data/Axis.h"
#include "hydrate/PlacementStrategy.h"
#include "hydrate/CullingStrategy.h"
#include "hydrate/GridMember.h"
#include "settings.h"
#include "Exceptions.h"

// forwards declaration
class Body;

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Grid {
  public:
    /**
     * @brief Construct a new Grid object with standard atomic radii.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     */
    Grid(const Axis3D& axes, double width) : Grid(axes, width, setting::grid::ra, setting::grid::rh, setting::grid::psc, setting::grid::csc) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     * @param radius the radius of each atom.
     */
    Grid(const Axis3D& axes, double width, int radius) : Grid(axes, width, radius, radius, setting::grid::psc, setting::grid::csc) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in each dimension. 
     * @param ra the radius of each atom.
     * @param rh the radius of each water molecule.
     */
    Grid(const Axis3D& axes, double width, double ra, double rh, setting::grid::PlacementStrategyChoice psc, setting::grid::CullingStrategyChoice csc);

    Grid(const vector<Atom>& atoms) : Grid(atoms, setting::grid::width, setting::grid::ra, setting::grid::rh, setting::grid::psc, setting::grid::csc) {}

    Grid(const vector<Atom>& atoms, double width, double ra, double rh, setting::grid::PlacementStrategyChoice psc, setting::grid::CullingStrategyChoice csc);

    /**
     * @brief Copy constructor. 
     */
    Grid(const Grid& grid);

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
     * @brief Add the contents of a body to this grid.
     * 
     * @param body The body to be added. 
     */
    vector<GridMember<Atom>> add(const Body* const body);

    /** 
     * @brief Add a single atom to the grid. 
     *        This is a constant-time operation. 
     */
    GridMember<Atom> add(const Atom& atom, const bool expand = false);

    /** 
     * @brief Add a single hetatom to the grid. 
     *        This is a constant-time operation. 
     */
    GridMember<Hetatom> add(const Hetatom& atom, const bool expand = false);

    /**
     * @brief Remove the contents of a body from this grid. 
     * 
     * @param body The body to be removed. 
     */
    void remove(const Body* const body);

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
     *        This is linear in the number of stored elements. 
     */
    void remove(const vector<Hetatom>& atom);

    /**
     * @brief Remove multiple protein atoms from the grid.
     *        This is linear in the number of stored elements. 
     */
    void remove(const vector<Atom>& atom);

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

    static std::pair<Vector3, Vector3> bounding_box(const vector<Atom>& atoms);

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

    vector<int> get_bins() const;

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
     * @brief Get the width of each bin.
     */
    double get_width() const {return width;}

    Axis3D get_axes() const {return axes;}

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

    /**
     * @brief Get a copy of this Grid. 
     */
    Grid copy() const;

    /**
     * @brief Set this Grid equal to another.
     * 
     * @param rhs The new Grid. 
     */
    Grid& operator=(const Grid& rhs);

    /**
     * @brief Check if this Grid is identical to another. 
     *        Primarily intended to be used with tests.
     * 
     * @param rhs The Grid to compare against. 
     */
    bool operator==(const Grid& rhs) const;

    vector<vector<vector<char>>> grid; // The actual grid. Datatype is char since we need at least four different values.
    std::list<GridMember<Atom>> a_members; // A list of all member atoms and where they are located.
    std::list<GridMember<Hetatom>> w_members; // A list of all member water molecules and where they are located. 
    int volume = 0; // The number of bins covered by the members, i.e. the actual volume in the unit (width)^3
    int ra = 0; // Radius of each atom represented as a number of bins
    int rh = 0; // Radius of each water molecule represented as a number of bins

  private:
    Axis3D axes;
    double width; // distance between each grid point
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

    void setup(double width, double ra, double rh, setting::grid::PlacementStrategyChoice psc, setting::grid::CullingStrategyChoice csc);
};