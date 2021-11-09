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
    Grid(TVector3 base, double width, int bins) : Grid(base, width, {bins, bins, bins}) {}

    Grid(TVector3 base, double width, vector<int> bins) {
        this->base = base;
        this->width = width;
        this->bins = bins;
        this->grid = vector(bins[0], vector<vector<bool>>(bins[1], vector<bool>(bins[2], false)));
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
     * @brief Expand the member atoms into actual spheres of a given radius. 
     * @param radius the radius of the spherical expansion in Ångström. 
     */
    void expand_volume(double radius) {
        set_radius(radius);

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

        // create a new HOH molecule based on an index vector ### NO SERIAL ###
        auto new_HOH = [&] (vector<int> v) {
            double x = base.X() + v[0]*width;
            double y = base.Y() + v[1]*width;
            double z = base.Z() + v[2]*width;
            return std::make_shared<Atom>(Atom({x, y, z}, 1, "O", "HOH", members.size()+1));
        };

        int c = 0; // counter
        for (vector<int> v : hydration_slots) {
            c++;
            shared_ptr<Atom> a = new_HOH(v);
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
    vector<vector<int>> find_free_locs() const {
        // a quick check to verify there are no water molecules already present
        for (const auto& pair : members) {
            if (pair.first.get_name() == "HOH") {
                print_err("WARNING: Attempting to hydrate a grid which already contains water!");
            }
        }
        if (r == 0) {
            print_err("WARNING: Attempting to hydrate a grid of atoms without any volume! This is most likely an error.");
        }

        vector<vector<int>> bounds = bounding_box();
        // add 2*radius extra space in each direction
        for (int i = 0; i < 3; i++) {
            bounds[i][0] = std::max(bounds[i][0] - 2*r, 0);
            bounds[i][1] = std::min(bounds[i][1] + 2*r, bins[i]);
        }

        // loop over the minimum bounding box as found above
        vector<vector<int>> available_locs;
        for (int i = bounds[0][0]; i < bounds[0][1]; i+=r) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j+=r) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k+=r) {
                    // if this spot is part of an atom
                    if (grid[i][j][k]) {
                        // check x ± r
                        if (!grid[std::max(i-r, 0)][j][k]) available_locs.push_back({std::max(i-r, 0), j, k});
                        if (!grid[std::min(i+r, bins[0])][j][k]) available_locs.push_back({std::min(i+r, bins[0]), j, k});

                        // check y ± r
                        if (!grid[i][std::max(j-r, 0)][k]) available_locs.push_back({i, std::max(j-r, 0), k});
                        if (!grid[i][std::min(j+r, bins[1])][k]) available_locs.push_back({i, std::min(j+r, bins[1]), k});

                        // check z ± r
                        if (!grid[i][j][std::max(k-r, 0)]) available_locs.push_back({i, j, std::max(k-r, 0)});
                        if (!grid[i][j][std::min(k+r, bins[2])]) available_locs.push_back({i, j, std::min(k+r, bins[2])});
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
            print_err("ERROR: Calculating a boundary box for a grid with no members!");
            exit(1);
        }

        // initialize the bounds with values that will always be replaced in the following loop
        vector<vector<int>> box = {{bins[0]+1, -1}, {bins[1]+1, -1}, {bins[2]+1, -1}};
        for (const auto& pair : members) {
            vector<int> loc = pair.second;
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
     * @brief Set the radius.
     * @param radius The new radius.
     */
    void set_radius(double radius) {
        int new_r = int(radius/width);
        if (this->r != 0 && this->r != new_r) {
            print_err("WARNING: The radius is already set for this grid!");
        }
        this->r = new_r;
    }

    /**
     * @brief Get all hydration atoms from this grid. 
     * NOTE: Consider adding a parameter to start all their serials from. 
     * @return A vector containing all of the hydration atoms. 
     */
    vector<Atom*> get_hydration_atoms() {
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
    int r = 0; // radius of each atom represented as a number of bins

    /** 
     * @brief Expand a single member atom into an actual sphere. The radius r must be set.
     */
    void expand_volume(Atom atom) {
        vector<int> loc = members.at(atom);
        for (int i = std::max(loc[0] - r, 0); i < std::min(loc[0] + r, bins[0]); i++) {
            for (int j = std::max(loc[1] - r, 0); j < std::min(loc[1] + r, bins[1]); j++) {
                for (int k = std::max(loc[2] - r, 0); k < std::min(loc[2] + r, bins[2]); k++) {
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
        int binx = (atom->get_x() - base.X()/width);
        int biny = (atom->get_y() - base.Y()/width);
        int binz = (atom->get_z() - base.Z()/width);
        members.insert({*atom, {binx, biny, binz}});
        grid[binx][biny][binz] = true;
    }
};