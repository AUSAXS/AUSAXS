#include <em/Datatypes.h>
#include <utility/Exceptions.h>
#include <ImageStack.h>
#include <settings/All.h>
#include <utility/Utility.h>
#include <crystal/Box.h>

#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
    settings::em::sample_frequency = 1;
    std::string file = argv[1];
    std::string output = "figures/" + utility::stem(file) + "/";
    utility::create_directory(output);

    em::ImageStack images(file);
    Box box = images.get_box(0);
    box.save(output + "box.dat");
    return 0;
}