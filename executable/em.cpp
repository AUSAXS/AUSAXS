#include <iostream>

#include <em/image.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    em::ImageStack image("data/A2M_map.ccp4");
    image.plot(std::stoi(argv[1]));
    return 0;
}