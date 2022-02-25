#include <iostream>

#include <em/image.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    em::ImageStack image("data/A2M_map.ccp4");
    image.plot(std::stoi(argv[1]));
    image.fit("data/A2M_ma.RSR");
    return 0;
}