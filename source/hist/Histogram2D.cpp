#include <hist/Histogram2D.h>

using namespace hist;

std::string Histogram2D::to_string() const {
    std::string str;
    str += "x_axis: " + std::to_string(x_axis.bins) + " " + std::to_string(x_axis.min) + " " + std::to_string(x_axis.max) + "\n";
    str += "x_axis: " + std::to_string(y_axis.bins) + " " + std::to_string(y_axis.min) + " " + std::to_string(y_axis.max) + "\n";
    for (unsigned int i = 0; i < x_axis.bins; i++) {
        for (unsigned int j = 0; j < y_axis.bins; j++) {
            str += std::to_string((*this)[i][j]) + " ";
        }
        str += "\n";
    }
    return str;
}