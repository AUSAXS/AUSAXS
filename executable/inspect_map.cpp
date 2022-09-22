#include <em/Datatypes.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <iostream>

int main(int argc, char const *argv[]) {
    //* Read existing header
    std::fstream file(argv[1], std::ios::in | std::ios::out | std::ios::binary);
    if (!file.is_open()) {throw except::io_error("Error in ImageStack::ImageStack: Could not open file \"" + std::string(argv[1]) + "\"");}

    em::ccp4::Header header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    std::cout << header << std::endl;

    return 0;
}