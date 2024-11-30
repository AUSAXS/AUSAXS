/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <plots/Styles.h>

#include <array>
#include <sstream>
#include <iomanip>

using namespace ausaxs;

std::string style::color_map::Rainbow::next() {
    std::array<double, 3> start = {1, 1, 1};

    // implementation adapted from https://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
    double v = static_cast<double>(i)/n;
    if (v < 0.25) {
        start[0] = 0;
        start[1] = 4*v;
    } else if (v < 0.5) {
        start[0] = 0;
        start[2] = 1 + 4*(0.25 - v);
    } else if (v < 0.75) {
        start[0] = 4*(v - 0.5);
        start[2] = 0;
    } else {
        start[1] = 1 + 4*(0.75 - v);
        start[2] = 0;
    }

    if (++i == n) {i = 0;}

    std::stringstream ss;
    ss << "#" << std::hex << std::setfill('0') 
        << std::setw(2) << static_cast<int>(start[0]*255) 
        << std::setw(2) << static_cast<int>(start[1]*255) 
        << std::setw(2) << static_cast<int>(start[2]*255);
    return ss.str();
}

std::string style::color_map::RedBlue::next() {
    std::array<double, 3> start = {153, 0, 0}, end = {0, 0, 153};

    double v = static_cast<double>(i)/n;
    std::array<double, 3> rgb = {
        start[0]*(1-v) + end[0]*v, 
        start[1]*(1-v) + end[1]*v, 
        start[2]*(1-v) + end[2]*v
    };

    std::stringstream ss;
    ss << "#" << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(rgb[0]) << std::setw(2) << static_cast<int>(rgb[1]) << std::setw(2) << static_cast<int>(rgb[2]);
    return ss.str();
}