#include <em/Datatypes.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <iostream>

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        std::cout << "Invalid number of arguments. First argument must be the map file, and the second the orientation in the form \"312\"." << std::endl;
        exit(1);
    }

    auto validate = [] (std::string order) {
        bool a = false, b = false, c = false;
        for (char v : order) {
            switch (v) {
                case '1': a = true; break;
                case '2': b = true; break;
                case '3': c = true; break;
                default: throw except::invalid_argument("Invalid argument. Received \"" + std::string(1, c) + "\", expected one of {1, 2, 3}.");
            }
        }
        return a && b && c;
    };

    std::string order = argv[2];
    if (order.size() != 3 || !validate(order)) {
        std::cout << "Invalid orientation form. Must be a permutation of \"123\". Received \"" + order + "\"" << std::endl;
        exit(1);
    }

    //* Read existing header
    std::fstream file(argv[1], std::ios::in | std::ios::out | std::ios::binary);
    if (!file.is_open()) {throw except::io_error("Error in ImageStack::ImageStack: Could not open file \"" + std::string(argv[1]) + "\"");}

    em::ccp4::Header header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));

    //* Rotate data
    int x = order[0] - '0';
    int y = order[1] - '0';
    int z = order[2] - '0';
    std::cout << "Rotating map from {" << header.mapc << ", " << header.mapr << ", " << header.maps << "} to {" << x << ", " << y << ", " << z << "}." << std::endl;
    header.rotate(x, y, z);

    //* Overwrite existing header
    char* data = reinterpret_cast<char*>(&header);
    file.seekp(0);
    file.write(data, sizeof(header));
    file.close();

    return 0;
}