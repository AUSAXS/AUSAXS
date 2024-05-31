/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/RadialLineGenerator.h>
#include <grid/Grid.h>
#include <constants/ConstantsMath.h>
#include <constants/vdwTable.h>
#include <math/Vector3.h>

using namespace grid::detail;

RadialLineGenerator::RadialLineGenerator(observer_ptr<grid::Grid> grid) {
    if (grid->get_width() != width) {
        width = grid->get_width();
        generate();
    }
}

void RadialLineGenerator::generate(int divisions) {
    std::vector<Vector3<int>> bins_1rh;
    std::vector<Vector3<int>> bins_3rh;
    std::vector<Vector3<int>> bins_5rh;
    std::vector<Vector3<int>> bins_7rh;
    std::vector<Vector3<double>> locs;
    double ang = 2*constants::pi/divisions;

    // we generate one octant of a sphere, and then reflect it to generate the rest
    // we do this to ensure the sphere is symmetric. If we simply generate it all at once, floating-point errors moves some of the bins around
    std::vector<Vector3<double>> sphere;
    for (double theta = 0; theta <= constants::pi*0.5; theta+=ang) {
        for (double phi = 0; phi <= constants::pi*0.5; phi+=ang) {
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

    double rh = constants::radius::vdw::O;
    for (const auto& rot : rots) {
        double xr = rot.x(), yr = rot.y(), zr = rot.z();
        bins_1rh.push_back(Vector3<int>(std::round(  rh*xr)/width, std::round(  rh*yr)/width, std::round(  rh*zr)/width));
        bins_3rh.push_back(Vector3<int>(std::round(3*rh*xr)/width, std::round(3*rh*yr)/width, std::round(3*rh*zr)/width));
        bins_5rh.push_back(Vector3<int>(std::round(5*rh*xr)/width, std::round(5*rh*yr)/width, std::round(5*rh*zr)/width));
        bins_7rh.push_back(Vector3<int>(std::round(7*rh*xr)/width, std::round(7*rh*yr)/width, std::round(7*rh*zr)/width));
        locs.push_back(Vector3<double>(xr, yr, zr));
    }

    // set the member vectors
    rot_bins_1rh = std::move(bins_1rh);
    rot_bins_3rh = std::move(bins_3rh);
    rot_bins_5rh = std::move(bins_5rh);
    rot_bins_7rh = std::move(bins_7rh);
    rot_locs = std::move(locs);
}