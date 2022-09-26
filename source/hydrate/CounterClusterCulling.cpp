#include <hydrate/CounterClusterCulling.h>
#include <hydrate/Grid.h>

using std::vector;

grid::CounterClusterCulling::CounterClusterCulling(Grid* grid) : ClusterCullingStrategy(grid) {
    prepare_rotations();
}

void grid::CounterClusterCulling::prepare_rotations() {
    unsigned int divisions = 8;

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

    int ra = grid->ra;
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

    std::unordered_map<unsigned int, unsigned int> groups; // maps location id to group id
    std::vector<std::vector<unsigned int>> group_members; // list of members in each group

    // reassign all members of group old_id to group new_id
    auto merge_groups = [&groups, &group_members] (unsigned int old_id, unsigned int new_id) {
        for (unsigned int member : group_members[old_id]) {
            groups[member] = new_id;
            group_members[new_id].push_back(member);
        }
        group_members[old_id].clear();
    };

    // add a member to a group
    auto add_to_group = [&groups, &group_members] (unsigned int id, unsigned int group_id) {
        groups[id] = group_id;
        group_members[group_id].push_back(id);
    };

    // create a new group with a single member
    auto add_to_new_group = [&groups, &group_members] (unsigned int id) {
        unsigned int group_id = group_members.size();
        groups[id] = group_id;
        group_members.push_back(std::vector<unsigned int>());
        group_members[group_id].push_back(id);
    };

    // Iterate through all atoms, and use radial lines to detect other nearby atoms. Group them into clusters.
    for (grid::GridMember<Atom>& atom : grid->a_members) {
        // check if atom is already in a group
        unsigned int id1 = to_id(atom.loc);

        // each atom starts in a group of its own, unless already added by someone else
        if (groups.count(id1) == 0) {
            add_to_new_group(id1);
        }

        // check spherical shell within 2ra for collisions
        for (const auto& bin : rot_bins_2ra) {
            Vector3<int> pos = atom.loc + bin;
            if (grid->grid.index(pos) == GridObj::EMPTY) {
                continue;
            }
            Vector3<int> center = find_center(pos);
            unsigned int id2 = to_id(center);

            // if the other atom is not in a group, add it to this one
            if (groups.count(id2) == 0) {
                add_to_group(id2, groups[id1]);
            } 
            
            // otherwise merge their group into ours
            else if (groups[id2] != groups[id1]) {
                merge_groups(groups[id2], groups[id1]);
            }
        }
    }

    // determine which groups to remove
    unsigned int sum = 0;
    unsigned int remove_count = 0;
    std::vector<bool> groups_to_remove(group_members.size(), false);
    for (unsigned int i = 0; i < group_members.size(); i++) {
        if (group_members[i].empty()) {continue;}
        if (group_members[i].size() < min_group_size) {
            groups_to_remove[i] = true;
            remove_count += group_members[i].size();
        }
        sum += group_members[i].size();
    }

    // sanity check
    if (sum != grid->a_members.size()) {
        throw except::unexpected("Error in CounterClusterCulling::cull: Group sizes (" + std::to_string(sum) + ") do not add up to total number of atoms (" + std::to_string(grid->a_members.size()) + ").");
    }

    // mark atoms for removal
    std::cout << remove_count << " atoms will be removed. " << std::endl;
    std::vector<bool> atoms_to_remove(grid->a_members.size(), false);
    unsigned int i = 0;
    for (auto& atom : grid->a_members) {
        unsigned int id = to_id(atom.loc);
        if (groups_to_remove[groups.at(id)]) {
            atoms_to_remove[i] = true;
            remove_count--;
        }
        i++;

        // if (atom.loc.x() > 1000 || atom.loc.y() > 1000 || atom.loc.z() > 1000) {
        //     std::cout << "Suspicious atom at " << atom.loc.to_string() << std::endl;
        // }
    }

    // sanity check
    if (remove_count != 0) {
        throw except::unexpected("Error in CounterClusterCulling::cull: Could not find all " + std::to_string(remove_count) + " atoms to be removed.");
    }

    return atoms_to_remove;
}