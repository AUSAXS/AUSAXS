#include <hydrate/CounterClusterCulling.h>
#include <hydrate/Grid.h>

using std::vector;

grid::CounterClusterCulling::CounterClusterCulling(Grid* grid) : ClusterCullingStrategy(grid) {
    prepare_rotations();
}

void grid::CounterClusterCulling::prepare_rotations() {
    unsigned int divisions = 6;

    int rh = grid->rh, ra = grid->ra;
    double width = grid->get_width();

    std::vector<Vector3<int>> bins_2ra;
    std::vector<Vector3<int>> bins_3ra;

    double ang = 2*M_PI/divisions;

    // we generate one octant of a sphere, and then reflect it to generate the rest
    // we do this to ensure the sphere is symmetric. If we simply generate it all at once, floating-point errors moves some of the bins around
    vector<Vector3<double>> sphere;
    for (double theta = 0; theta <= M_PI*0.5; theta+=ang) {
        for (double phi = 0; phi <= M_PI*0.5; phi+=ang) {
            double x = cos(phi)*sin(theta);
            double y = sin(phi)*sin(theta);
            double z = cos(theta);
            sphere.push_back({x, y, z});
            sphere.push_back({-x, y, z});
            sphere.push_back({x, -y, z});
            sphere.push_back({-x, -y, z});
            sphere.push_back({x, y, -z});
            sphere.push_back({-x, y, -z});
            sphere.push_back({x, -y, -z});
            sphere.push_back({-x, -y, -z});
        }
    }

    // remove duplicates
    vector<Vector3<double>> rots;
    for (auto& p : sphere) {
        bool present = false;
        for (int i = 0; i < 3; i++) { // fix the easy floating point errors
            if (abs(p[i]) < 1e-5) {p[i] = 0;}
        }
        for (const auto& r : rots) { // go through all rotations and try to find a duplicate entry
            if (r.distance(p) < 1e-5) {
                present = true;
                break;
            }
        }
        if (!present) { // if the element was not already present
            rots.push_back(p); // add it
        }
    }

    double rarh = ra+rh;
    for (const auto& rot : rots) {
        double xr = rot.x(), yr = rot.y(), zr = rot.z();
        bins_2ra.push_back(Vector3<int>(std::trunc(2*ra*xr), std::trunc(2*ra*yr), std::trunc(2*ra*zr)));
        bins_3ra.push_back(Vector3<int>(std::trunc(3*ra*xr), std::trunc(3*ra*yr), std::trunc(3*ra*zr)));
    }

    // set the member vectors
    rot_bins_2ra = std::move(bins_2ra);
    rot_bins_3ra = std::move(bins_3ra);
}

vector<bool> grid::CounterClusterCulling::cull(unsigned int min_group_size) const {
    // find the center of an expanded atom
    unsigned int ra = grid->ra;
    auto find_center = [this, ra] (const Vector3<int> pos) {
        // we split the search for each axis into two parts so we can break out of the loop early if we leave the atom
        // so instead of (pos_x - ra) to (pos_x + ra), we do (pos_x) to (pos_x - ra) and (pos_x) to (pos_x + ra) (repeat for y & z)

        // [xmin, x]
        for (unsigned int x = pos.x(); x >= pos.x()-ra; x--) {
            // [ymin, y]
            for (unsigned int y = pos.y(); y >= pos.y()-ra; y--) {
                // [zmin, z]
                for (unsigned int z = pos.z(); z >= pos.z()-ra; z--) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }

                // [z, zmax]
                for (unsigned int z = pos.z(); z <= pos.z()+ra; z++) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }
            }

            // [y, ymax]
            for (unsigned int y = pos.y(); y <= pos.y()+ra; y++) {
                // [zmin, z]
                for (unsigned int z = pos.z(); z >= pos.z()-ra; z--) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }

                // [z, zmax]
                for (unsigned int z = pos.z(); z <= pos.z()+ra; z++) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }
            }
        }

        // [x, xmax]
        for (unsigned int x = pos.x(); x <= pos.x()+ra; x++) {
            // [ymin, y]
            for (unsigned int y = pos.y(); y >= pos.y()-ra; y--) {
                // [zmin, z]
                for (unsigned int z = pos.z(); z >= pos.z()-ra; z--) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }

                // [z, zmax]
                for (unsigned int z = pos.z(); z <= pos.z()+ra; z++) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }
            }

            // [y, ymax]
            for (unsigned int y = pos.y(); y <= pos.y()+ra; y++) {
                // [zmin, z]
                for (unsigned int z = pos.z(); z >= pos.z()-ra; z--) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }

                // [z, zmax]
                for (unsigned int z = pos.z(); z <= pos.z()+ra; z++) {
                    auto cell = grid->grid.index(x, y, z);
                    if (cell == GridObj::EMPTY) {break;}
                    else if (cell == GridObj::A_CENTER) {return Vector3<int>(x, y, z);}
                }
            }
        }

        throw except::unexpected("Could not find center of atom at " + pos.to_string() + " in grid.");
    };

    // convert a position to a unique id
    auto to_id = [] (const Vector3<int> pos) {
        return pos.x()*1e8 + pos.y()*1e4 + pos.z();
    };

    // Iterate through all atoms, and use radial lines to detect other nearby atoms. Group them into clusters.
    std::unordered_map<unsigned int, unsigned int> groups; // maps location id to group id
    std::vector<unsigned int> group_sizes; // number of members in each group
    for (grid::GridMember<Atom>& atom : grid->a_members) {
        // check if atom is already in a group
        unsigned int id1 = to_id(atom.loc);
        if (groups.count(id1) != 0) {
            continue;
        }

        // check spherical shell within 2ra
        for (const auto& bin : rot_bins_2ra) {
            Vector3<int> pos = atom.loc + bin;
            if (grid->grid.index(pos) == GridObj::EMPTY) {
                continue;
            }
            Vector3<int> center = find_center(pos);
            unsigned int id2 = to_id(center);

            // check if the other atom is already in a group
            if (groups.count(id2) == 0) {
                // if not, create a new one for them
                unsigned int group_no = group_sizes.size();
                groups[id1] = group_no;
                groups[id2] = group_no;
                group_sizes.push_back(2);
            } 
            
            else {
                // otherwise add this atom to the same group
                groups[id1] = groups.at(id2);
                group_sizes[groups.at(id2)]++;
            }
        }
    }

    // determine which groups to remove
    std::vector<bool> groups_to_remove(group_sizes.size(), false);
    for (unsigned int i = 0; i < group_sizes.size(); i++) {
        if (group_sizes[i] < min_group_size) {
            groups_to_remove[i] = true;
        }
    }

    // mark atoms for removal
    std::vector<bool> atoms_to_remove(grid->a_members.size(), false);
    unsigned int i = 0;
    for (auto& atom : grid->a_members) {
        unsigned int id = to_id(atom.loc);
        if (groups_to_remove[groups.at(id)]) {
            atoms_to_remove[i] = true;
        }
        i++;
    }

    return atoms_to_remove;
}