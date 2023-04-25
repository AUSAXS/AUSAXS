#include <em/Box.h>
#include <utility/Exceptions.h>
#include <utility/Axis.h>

#include <fstream>
#include <numeric>

using namespace crystal;

void Box::save(std::string filename) const {
    std::ofstream output(filename);
    if (!output.is_open()) {throw except::io_error("Box::save: Could not open file \"" + filename + "\"");}

    output << "x y z" << std::endl;
    for (const auto& atom : atoms) {
        output << atom.x() << " " << atom.y() << " " << atom.z() << std::endl;
    }
}

CoordinateSystem Box::spanning_coordinate_system() const {
    // find the eight corner atoms
    Vector3<float> min_x_min_y_min_z, min_x_min_y_max_z, min_x_max_y_min_z, min_x_max_y_max_z, max_x_min_y_min_z, max_x_min_y_max_z, max_x_max_y_min_z, max_x_max_y_max_z;
    for (const auto& atom : atoms) {
        if (atom.x() < min_x_min_y_min_z.x() && atom.y() < min_x_min_y_min_z.y() && atom.z() < min_x_min_y_min_z.z()) {min_x_min_y_min_z = atom;}
        if (atom.x() < min_x_min_y_max_z.x() && atom.y() < min_x_min_y_max_z.y() && atom.z() > min_x_min_y_max_z.z()) {min_x_min_y_max_z = atom;}
        if (atom.x() < min_x_max_y_min_z.x() && atom.y() > min_x_max_y_min_z.y() && atom.z() < min_x_max_y_min_z.z()) {min_x_max_y_min_z = atom;}
        if (atom.x() < min_x_max_y_max_z.x() && atom.y() > min_x_max_y_max_z.y() && atom.z() > min_x_max_y_max_z.z()) {min_x_max_y_max_z = atom;}
        if (atom.x() > max_x_min_y_min_z.x() && atom.y() < max_x_min_y_min_z.y() && atom.z() < max_x_min_y_min_z.z()) {max_x_min_y_min_z = atom;}
        if (atom.x() > max_x_min_y_max_z.x() && atom.y() < max_x_min_y_max_z.y() && atom.z() > max_x_min_y_max_z.z()) {max_x_min_y_max_z = atom;}
        if (atom.x() > max_x_max_y_min_z.x() && atom.y() > max_x_max_y_min_z.y() && atom.z() < max_x_max_y_min_z.z()) {max_x_max_y_min_z = atom;}
        if (atom.x() > max_x_max_y_max_z.x() && atom.y() > max_x_max_y_max_z.y() && atom.z() > max_x_max_y_max_z.z()) {max_x_max_y_max_z = atom;}
    }

    // find the two corner atoms with the smallest distance between them
    double min_distance = std::numeric_limits<double>::max();
    Vector3<float> corner1, corner2;
    for (const auto& atom1 : {min_x_min_y_min_z, min_x_min_y_max_z, min_x_max_y_min_z, min_x_max_y_max_z, max_x_min_y_min_z, max_x_min_y_max_z, max_x_max_y_min_z, max_x_max_y_max_z}) {
        for (const auto& atom2 : {min_x_min_y_min_z, min_x_min_y_max_z, min_x_max_y_min_z, min_x_max_y_max_z, max_x_min_y_min_z, max_x_min_y_max_z, max_x_max_y_min_z, max_x_max_y_max_z}) {
            if (atom1 == atom2) {continue;}
            double distance = (atom1 - atom2).norm();
            if (distance < min_distance) {
                min_distance = distance;
                corner1 = atom1;
                corner2 = atom2;
            }
        }
    }

    Vector3 axis1 = corner2 - corner1;

    // define the second axis as the one with the largest angle to the first
    Vector3<float> axis2, axis3;
    double max_angle = 0;
    for (const auto& atom : atoms) {
        Vector3 axis = atom - corner1;
        double angle = std::acos(axis1.dot(axis) / (axis1.norm() * axis.norm()));
        if (angle > max_angle) {
            max_angle = angle;
            axis2 = axis;
        }
    }

    // define the third axis as the one with the largest angle to the two others
    max_angle = 0;
    for (const auto& atom : atoms) {
        Vector3<float> axis = atom - corner1;
        double angle = std::acos(axis1.dot(axis) / (axis1.norm() * axis.norm())) + std::acos(axis2.dot(axis) / (axis2.norm() * axis.norm()));
        if (angle > max_angle) {
            max_angle = angle;
            axis3 = axis;
        }
    }

    // normalize axes
    // axis1.normalize();
    // axis2.normalize();
    // axis3.normalize();

    return CoordinateSystem(corner1, axis1, axis2, axis3);
}

void Box::distances() const {
    // generous sizes - 2000Å should be enough for just about any structure
    double width = 1;
    Axis axes = Axis(0, 5000, 5000/width); 
    std::vector<double> p_pp(axes.bins, 0);

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time
    std::vector<float> data_p(atoms.size()*4);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        const auto& a = atoms[i]; 
        data_p[4*i] = a.x();
        data_p[4*i+1] = a.y();
        data_p[4*i+2] = a.z();
        data_p[4*i+3] = 1;
    }

    // calculate p-p distances
    for (unsigned int i = 0; i < atoms.size(); i++) {
        for (unsigned int j = i+1; j < atoms.size(); j++) {
            float weight = data_p[4*i+3]*data_p[4*j+3]; // Z1*Z2*w1*w2
            float dx = data_p[4*i] - data_p[4*j];
            float dy = data_p[4*i+1] - data_p[4*j+1];
            float dz = data_p[4*i+2] - data_p[4*j+2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp.at(dist/width) += 2*weight;
            // p_pp[dist/width] += 2*weight;
        }
    }

    // add self-correlation
    for (unsigned int i = 0; i < atoms.size(); i++) {
        p_pp[0] += data_p[4*i+3]*data_p[4*i+3];
    }

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
    for (int i = axes.bins-1; i >= 10; i--) {
        if (p_pp[i] != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_pp.resize(max_bin);
    axes = Axis{max_bin, 0, max_bin*width}; 

    // write to file
    std::ofstream file;
    file.open("distances.dat");
    for (int i = 0; i < axes.bins; i++) {
        file << i << " " << p_pp[i] << std::endl;
    }
}