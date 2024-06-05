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

std::vector<Vector3<int>> RadialLineGenerator::rot_bins_1;
std::vector<Vector3<int>> RadialLineGenerator::rot_bins_2;
std::vector<Vector3<int>> RadialLineGenerator::rot_bins_3;
std::vector<Vector3<int>> RadialLineGenerator::rot_bins_4;
std::vector<Vector3<double>> RadialLineGenerator::rot_locs_abs;

RadialLineGenerator::RadialLineGenerator(observer_ptr<grid::Grid> grid, double radius, int divisions) : RadialLineGenerator(grid, {radius, 3*radius, 5*radius, 7*radius}, divisions) {}

RadialLineGenerator::RadialLineGenerator(observer_ptr<grid::Grid> grid, std::array<double, 4> radius, int divisions) {
    if (grid->get_width() != width) {
        width = grid->get_width();
        generate(radius, divisions);
    }
}

RadialLineGenerator::~RadialLineGenerator() = default;

void RadialLineGenerator::generate(std::array<double, 4> r, int divisions) {
    std::vector<Vector3<int>> bins_1, bins_2, bins_3, bins_4;
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

    for (const auto& rot : rots) {
        double xr = rot.x(), yr = rot.y(), zr = rot.z();
        bins_1.push_back(Vector3<int>(std::round(r[0]*xr)/width, std::round(r[0]*yr)/width, std::round(r[0]*zr)/width));
        bins_2.push_back(Vector3<int>(std::round(r[1]*xr)/width, std::round(r[1]*yr)/width, std::round(r[1]*zr)/width));
        bins_3.push_back(Vector3<int>(std::round(r[2]*xr)/width, std::round(r[2]*yr)/width, std::round(r[2]*zr)/width));
        bins_4.push_back(Vector3<int>(std::round(r[3]*xr)/width, std::round(r[3]*yr)/width, std::round(r[3]*zr)/width));
        locs.push_back(Vector3<double>(xr, yr, zr));
    }

    // set the member vectors
    rot_bins_1 = std::move(bins_1);
    rot_bins_2 = std::move(bins_2);
    rot_bins_3 = std::move(bins_3);
    rot_bins_4 = std::move(bins_4);
    rot_locs_abs = std::move(locs);
}