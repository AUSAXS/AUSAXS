#include <Symbols.h>
#include <hydrate/RadialPlacement.h>
#include <hydrate/Grid.h>
#include <utility/Settings.h>

void grid::RadialPlacement::prepare_rotations(int divisions) {
    int rh = grid->get_radius_water(), ra = grid->get_radius_atoms();
    double width = grid->get_width();

    // vector<vector<int>> bins;
    std::vector<Vector3<int>> bins_1rh;
    std::vector<Vector3<int>> bins_3rh;
    std::vector<Vector3<int>> bins_5rh;
    std::vector<Vector3<int>> bins_7rh;
    std::vector<Vector3<int>> bins_rarh;
    std::vector<Vector3<double>> locs_rarh;
    double ang = 2*M_PI/divisions;

    // we generate one octant of a sphere, and then reflect it to generate the rest
    // we do this to ensure the sphere is symmetric. If we simply generate it all at once, floating-point errors moves some of the bins around
    std::vector<Vector3<double>> sphere;
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
    std::vector<Vector3<double>> rots;
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
        bins_1rh.push_back(Vector3<int>(std::trunc(rh*xr), std::trunc(rh*yr), std::trunc(rh*zr)));
        bins_3rh.push_back(Vector3<int>(std::trunc(3*rh*xr), std::trunc(3*rh*yr), std::trunc(3*rh*zr)));
        bins_5rh.push_back(Vector3<int>(std::trunc(5*rh*xr), std::trunc(5*rh*yr), std::trunc(5*rh*zr)));
        bins_7rh.push_back(Vector3<int>(std::trunc(7*rh*xr), std::trunc(7*rh*yr), std::trunc(7*rh*zr)));
        bins_rarh.push_back(Vector3<int>(std::trunc((rarh)*xr), std::trunc((rarh)*yr), std::trunc((rarh)*zr)));
        locs_rarh.push_back(Vector3<double>((rarh*width)*xr, (rarh*width)*yr, (rarh*width)*zr));
    }

    // set the member vectors
    rot_bins_1rh = std::move(bins_1rh);
    rot_bins_3rh = std::move(bins_3rh);
    rot_bins_5rh = std::move(bins_5rh);
    rot_bins_7rh = std::move(bins_7rh);
    rot_bins_rarh = std::move(bins_rarh);
    rot_locs_rarh = std::move(locs_rarh);
}

std::vector<grid::GridMember<Water>> grid::RadialPlacement::place() const {
    // dereference the values we'll need for better performance
    auto bins = grid->get_bins();
    GridObj& gref = grid->grid;

    // we define a helper lambda
    std::vector<GridMember<Water>> placed_water(grid->a_members.size());
    unsigned int index = 0;
    auto add_loc = [&] (Vector3<double> exact_loc) {
        Water a = Water::create_new_water(exact_loc);
        GridMember<Water> gm = grid->add(a, true);
        if (placed_water.size() <= index) [[unlikely]] {
            placed_water.resize(2*index);
        }
        placed_water[index++] = gm;
    };

    std::cout << "Atom order: " << std::endl;
    for (const auto& atom : grid->a_members) {
        std::cout << atom.loc << std::endl;
    }

    for (const auto& atom : grid->a_members) {
        int x = atom.loc.x(), y = atom.loc.y(), z = atom.loc.z();
        std::cout << "Placing water around atom at location " << x << " " << y << " " << z << std::endl;

        for (unsigned int i = 0; i < rot_bins_rarh.size(); i++) {
            int xr = x + rot_bins_rarh[i].x(), yr = y + rot_bins_rarh[i].y(), zr = z + rot_bins_rarh[i].z(); // new coordinates
            
            // check bounds
            if (xr < 0) xr = 0;
            if (xr >= (int) bins.x()) xr = bins.x()-1;
            if (yr < 0) yr = 0;
            if (yr >= (int) bins.y()) yr = bins.y()-1;
            if (zr < 0) zr = 0;
            if (zr >= (int) bins.z()) zr = bins.z()-1;

            std::cout << "\tChecking location " << xr << " " << yr << " " << zr << std::endl;
            // we have to make sure we don't check the direction of the atom we are trying to place this water on
            Vector3<int> skip_bin(xr-rot_bins_1rh[i].x(), yr-rot_bins_1rh[i].y(), zr-rot_bins_1rh[i].z());
            if (gref.index(xr, yr, zr) == GridObj::EMPTY && collision_check(Vector3<int>(xr, yr, zr), skip_bin)) {
                Vector3<double> exact_loc = atom.atom.coords + rot_locs_rarh[i];
                add_loc(exact_loc);
            } else {
                std::cout << "\t\tInvalid location" << std::endl;
            }
        }
    }

    placed_water.resize(index);
    return placed_water;
}

bool grid::RadialPlacement::collision_check(const Vector3<int>& loc, const Vector3<int>& skip_bin) const {
    std::cout << "\t\tChecking location " << loc << std::endl;
    // dereference the values we'll need for better performance
    GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    int score = 0;
    // check if a location is out-of-bounds
    auto is_out_of_bounds = [&bins, &score] (Vector3<int> v) {
        if (v.x() < 0 || (int) bins.x() <= v.x() ) {return true;}
        if (v.y() < 0 || (int) bins.y() <= v.y() ) {return true;}
        if (v.z() < 0 || (int) bins.z() <= v.z() ) {return true;}
        return false;
    };

    for (unsigned int i = 0; i < rot_bins_1rh.size(); i++) {
        // check for collisions at 1rh
        int xr = loc.x() + rot_bins_1rh[i].x(), yr = loc.y() + rot_bins_1rh[i].y(), zr = loc.z() + rot_bins_1rh[i].z(); // new coordinates

        // check for bounds
        if (xr < 0) xr = 0;
        if (xr >= (int) bins.x()) xr = bins.x()-1;
        if (yr < 0) yr = 0;
        if (yr >= (int) bins.y()) yr = bins.y()-1;
        if (zr < 0) zr = 0;
        if (zr >= (int) bins.z()) zr = bins.z()-1;

        if (gref.index(xr, yr, zr) != GridObj::EMPTY) {
            if (Vector3(xr, yr, zr) == skip_bin) {continue;} // skip the bin containing the atom we're trying to place this water molecule on
            std::cout << "\t\t\tCollision at 1rh" << std::endl;
            return false;
        }

        // check if we're in a cavity
        for (unsigned int j = 0; j < rot_bins_3rh.size(); j++) {
            // check at 3r
            int xr = loc.x() + rot_bins_3rh[j].x(), yr = loc.y() + rot_bins_3rh[j].y(), zr = loc.z() + rot_bins_3rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) { // if the line goes out of bounds, we know for sure it won't intersect anything
                score += 3;                       // so we add three points and move on to the next
                continue;
            }
            if (gref.index(xr, yr, zr) != GridObj::EMPTY) { // if the line intersects something at 3r, we don't check the other two points of the same line
                score -= 3;                                 // but immediately subtract 3 points and move on to the next
                continue;
            } else {
                score++;
            } 

            // check at 5r
            xr = loc.x() + rot_bins_5rh[j].x(); yr = loc.y() + rot_bins_5rh[j].y(); zr = loc.z() + rot_bins_5rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 2;
                continue;
            }
            if (gref.index(xr, yr, zr) != GridObj::EMPTY) {
                score -= 2;
                continue;
            } else {
                score++;
            } 

            // check at 7r
            xr = loc.x() + rot_bins_7rh[j].x(); yr = loc.y() + rot_bins_7rh[j].y(); zr = loc.z() + rot_bins_7rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 1;
                continue;
            }
            if (gref.index(xr, yr, zr) != GridObj::EMPTY) {
                score -= 1;
                continue;
            } else {
                score++;
            }
        }
    }
    if (score <= setting::grid::placement::min_score*rot_bins_1rh.size()) {
        std::cout << "\t\t\tScored too low" << std::endl;
        return false;
    }
    std::cout << "\t\tAccepted" << std::endl;
    return true;
}