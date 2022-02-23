#include <iostream>

#include <em/image.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    std::ifstream input("data/A2M_map.ccp4", std::ios::binary);
    em::ccp4::Header header;
    input.read(reinterpret_cast<char*>(&header), sizeof(header));

    em::Image image(header, input);
    image.plot(std::stoi(argv[1]));

    std::cout << header.nx << ", " << header.ny << ", " << header.nz << ", " << header.mode << std::endl;
    std::cout << header.cella_x << ", " << header.cella_y << ", " << header.cella_z << std::endl;
    std::cout << header.nlabl << std::endl;
    std::cout << string(header.label.data()).substr(0, 80) << std::endl;
    return 0;
}