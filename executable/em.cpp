#include <iostream>

#include <em/image.h>
#include <plots/PlotImage.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    // em::ImageStack image("data/A2M_map.ccp4");
    // image.plot(std::stoi(argv[1]));
    // image.fit("data/A2M_ma.RSR");

    em::ImageStack image("data/A2M_map.ccp4"); 
    PlotImage plot(image.image(5));
    plot.plot_atoms(0.2);
    plot.save("test.pdf");

    return 0;
}