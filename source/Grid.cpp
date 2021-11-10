#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "data/Atom.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class Grid {
public:
    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     */
    Grid(TVector3 base, double width, int bins) : Grid(base, width, {bins, bins, bins}, sqrt(8)) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in all dimensions. 
     * @param radius the radius of each atom.
     */
    Grid(TVector3 base, double width, int bins, int radius) : Grid(base, width, {bins, bins, bins}, radius) {}

    /**
     * @brief Construct a new Grid object.
     * @param base the base point for the grid.
     * @param width the distance between each point.
     * @param bins the number of bins in each dimension. 
     * @param radius the radius of each atom.
     */
    Grid(TVector3 base, double width, vector<int> bins, double radius) {
        this->base = base;
        this->width = width;
        this->bins = bins;
        this->grid = vector(bins[0], vector<vector<bool>>(bins[1], vector<bool>(bins[2], false)));
        this->set_radius_atoms(radius);
        this->set_radius_water(radius);
    }

    /** 
     * @brief Add a set of atoms to the grid. 
     * @param atoms the set of atoms to add to this grid.
     */
    void add(vector<shared_ptr<Atom>>* atoms) {
        for (auto const& a : *atoms) {
            add(a);
        }
    }

    /** 
     * @brief Expand the member atoms into actual spheres of radius r. 
     */
    void expand_volume() {
        vol_expanded = true;

        // iterate through each member location
        for (const auto& pair : members) {
            expand_volume(pair.first);
        }
    }

    /**
     * @brief Hydrate the grid. 
     * @param reduce Reduces the number of generated HOH molecules by this factor. Use 0 for no reduction. 
     */
    vector<shared_ptr<Atom>> hydrate(int reduce = 3) {
        vector<shared_ptr<Atom>> hydration_atoms;
        vector<vector<int>> hydration_slots = find_free_locs();

        // create a new water molecule based on an index vector ### NO SERIAL ###
        auto create_new_water = [&] (vector<int> v) {
            double x = base.X() + v[0]*width;
            double y = base.Y() + v[1]*width;
            double z = base.Z() + v[2]*width;
            return Atom::create_new_water({x, y, z});
        };

        int c = 0; // counter
        for (vector<int> v : hydration_slots) {
            c++;
            shared_ptr<Atom> a = create_new_water(v);
            if (reduce != 0) {
                if (c % reduce != 0) {
                    continue;
                }
            }
            add(a);
            expand_volume(*a);
            hydration_atoms.push_back(a);
        }

        return hydration_atoms;
    }

    /**
     * @brief Identify possible hydration binding locations for the structure. 
     * @return A list of possible (binx, biny, binz) locations.
     */
    vector<vector<int>> find_free_locs() {
        // a quick check to verify there are no water molecules already present
        for (const auto& pair : members) {
            if (pair.first.is_water()) {
                print_err("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
            }
        }

        if (!vol_expanded) {
            expand_volume();
        }

        vector<vector<int>> bounds = bounding_box();
        // add ra+rh extra space in each direction to account for the volume of both atoms and water molecules
        for (int i = 0; i < 3; i++) {
            bounds[i][0] = std::max(bounds[i][0] - (ra+rh), 0);
            bounds[i][1] = std::min(bounds[i][1] + (ra+rh) + 1, bins[i]); // +1 since this range is inclusive, but the following loop is not
        }

        // loop over the minimum bounding box as found above
        vector<vector<int>> available_locs;
        auto check_loc = [&] (const vector<int> v) {return check_collisions(v, available_locs);};
        for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                    // if this spot is part of an atom
                    if (grid[i][j][k]) {
                        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                        int xmin = std::max(i-rh, 0), xmax = std::min(i+rh, bins[0]);
                        int ymin = std::max(j-rh, 0), ymax = std::min(j+rh, bins[1]);
                        int zmin = std::max(k-rh, 0), zmax = std::min(k+rh, bins[2]);

                        // check collisions for x ± r_eff
                        if (!grid[xmin][j][k] && check_loc({xmin, j, k})) available_locs.push_back({xmin, j, k});
                        if (!grid[xmax][j][k] && check_loc({xmax, j, k})) available_locs.push_back({xmax, j, k});

                        // check collisions for y ± r_eff
                        if (!grid[i][ymin][k] && check_loc({i, ymin, k})) available_locs.push_back({i, ymin, k});
                        if (!grid[i][ymax][k] && check_loc({i, ymax, k})) available_locs.push_back({i, ymax, k});

                        // check collisions for z ± r_eff
                        if (!grid[i][j][zmin] && check_loc({i, j, zmin})) available_locs.push_back({i, j, zmin});
                        if (!grid[i][j][zmax] && check_loc({i, j, zmax})) available_locs.push_back({i, j, zmax});
                    }
                }
            }
        }
        cout << "Found " << available_locs.size() << " available HOH spots." << endl;
        return available_locs;
    };

    /**
     * @brief Create the smallest possible box containing the center points of all member atoms.
     * @return vector<vector<int>> An index pair (min, max) for each dimension (shape: [3][2]). 
     */
    vector<vector<int>> bounding_box() const {
        if (members.size() == 0) {
            print_err("Error in Grid::bounding_box: Calculating a boundary box for a grid with no members!");
            exit(1);
        }

        // initialize the bounds as large as possible
        vector<vector<int>> box = {{bins[0], 0}, {bins[1], 0}, {bins[2], 0}};
        for (const auto& pair : members) {
            const vector<int>& loc = pair.second;
            for (int i = 0; i < 3; i++) {
                if (box[i][0] > loc[i]) box[i][0] = loc[i]; // min
                if (box[i][1] < loc[i]) box[i][1] = loc[i]; // max
            }
        }
        return box;
    }

    /** 
     * @brief Returns a pointer to the boolean grid. 
     * @return A pointer to the boolean grid. 
     */
    vector<vector<vector<bool>>>* get_grid() {return &grid;}

    /**
     * @brief Set the radius of all atoms (not water molecules!).
     * @param radius The new radius in Ångström.
     */
    void set_radius_atoms(double radius) {
        int new_r = int(radius/width); // convert the radius to a "bin-radius"
        if (this->ra != 0 && this->ra != new_r) {
            print_err("Warning in Grid::set_radius: The radius is already set for this grid!");
        }
        this->ra = new_r;
    }

    /**
     * @brief Set the radius for all water molecules.
     * @param radius The new radius in Ångström.
     */
    void set_radius_water(double radius) {
        int new_r = int(radius/width); // convert the radius to a "bin-radius"
        if (this->rh != 0 && this->rh != new_r) {
            print_err("Warning in Grid::set_radius: The radius is already set for this grid!");
        }
        this->rh = new_r;
    }

    /**
     * @brief Get all hydration atoms from this grid. 
     * NOTE: Consider adding a parameter to start all their serials from. 
     * @return A vector containing all of the hydration atoms. 
     */
    vector<Atom*> get_hydration_atoms() const {
        vector<Atom*> atoms;
        for (const auto& pair : members) {
            Atom a = pair.first;
            if (a.is_water()) {
                atoms.push_back(&a);
            }
        }
        return atoms;
    }

private:
    TVector3 base; // base point of this grid
    double width; // distance between each grid point
    vector<vector<vector<bool>>> grid; // the actual grid
    std::map<Atom, vector<int>> members; // a map of all members of this grid and where they are located
    vector<int> bins; // the number of bins in each dimension
    int volume = 0; // the number of bins covered by the members, i.e. the actual volume in the unit (width)^3
    int ra = 0; // radius of each atom represented as a number of bins
    int rh = 0; // radius of each water molecule represented as a number of bins
    bool vol_expanded = false; // a flag determining if the volume has been expanded 

    /**
     * @brief Check if a water molecule can be placed at the given location. 
     *        This checks collisions with both other water molecules and other atoms. 
     * @param loc the location to be checked. 
     * @param other_molecules the other water molecules which have already been placed.
     * @return True if this is an acceptable location, false otherwise.
     */
    bool check_collisions(const vector<int> loc, const vector<vector<int>> other_molecules) const {
        const int x = loc[0], y = loc[1], z = loc[2];

        double amin = 100, wmin = 100;

        // check collision with other atoms
        int r_eff = ra + rh; // effective radius for water/atom collisions
        for (auto const& pair : members) {
            const vector<int>& v = pair.second;

            // check if the point is inside the spherical volume of loc
            double d = sqrt(pow(x-v[0], 2) + pow(y-v[1], 2) + pow(z-v[2], 2));
            amin = std::min(amin, d);
            if (d < r_eff) {
                // cout << format("Atom: (%1%, %2%, %3%), water: (%4%, %5%, %6%)") % v[0] % v[1] % v[2] % x % y % z << endl;
                // cout << "Separation: " << sqrt(pow(x-v[0], 2) + pow(y-v[1], 2) + pow(z-v[2], 2)) << endl << endl;
                return false;
            }
        }

        // check collision with other water molecules
        r_eff = 2*rh; // effective radius for water/water collisions
        for (auto const& a : other_molecules) {
            // check if the point is inside the spherical volume of loc
            double d = sqrt(pow(x-a[0], 2) + pow(y-a[1], 2) + pow(z-a[2], 2));
            wmin = std::min(wmin, d);
            if (sqrt(pow(x-a[0], 2) + pow(y-a[1], 2) + pow(z-a[2], 2)) < r_eff) {
                return false;
            }            
        }
        cout << "Water molecule placed! (amin, wmin): (" << amin << ", " << wmin << ")" << endl;
        return true;
    }

    /** 
     * @brief Expand a single member atom into an actual sphere.
     * @param atom the member atom to be expanded. 
     */
    void expand_volume(const Atom atom) {
        vector<int> loc = members.at(atom);
        int r = atom.is_water() ? rh : ra; // determine which radius to use for the expansion

        // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        vector<vector<int>> bounds(3, vector<int>(2, 0));
        for (int i = 0; i < 3; i++) {
            bounds[i][0] = std::max(loc[i] - r, 0);
            bounds[i][1] = std::min(loc[i] + r + 1, bins[i]); // +1 since this range is inclusive, while the following for-loop is not
        }

        // loop over each bin in the box
        for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                    // determine if the bin is within a sphere centered on the atom
                    if (std::sqrt(std::pow(loc[0] - i, 2) + std::pow(loc[1] - j, 2) + std::pow(loc[2] - k, 2)) <= r) {
                        volume++;
                        grid[i][j][k] = true;
                    }
                }
            }
        }
    }

    /** 
     * @brief Add a single atom to the grid. 
     * @param atoms The atom to be added. 
     */
    void add(shared_ptr<Atom> atom) {
        volume++;
        int binx = std::round(atom->get_x() - base.X()/width);
        int biny = std::round(atom->get_y() - base.Y()/width);
        int binz = std::round(atom->get_z() - base.Z()/width);
        members.insert({*atom, {binx, biny, binz}});
        grid[binx][biny][binz] = true;
    }
};