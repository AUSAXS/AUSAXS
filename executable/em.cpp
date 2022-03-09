#include <iostream>

#include <em/ImageStack.h>
#include <plots/PlotImage.h>
#include <plots/PlotIntensity.h>
#include <plots/PlotDistance.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    // em::ImageStack image("data/A2M_map.ccp4");
    // image.plot(std::stoi(argv[1]));
    // image.fit("data/A2M_ma.RSR");

    em::ImageStack image("data/A2M_map.ccp4"); 
    PlotImage plot(image.image(std::stoi(argv[1])));
    plot.plot_atoms(-1);
    plot.save("temp.pdf");
    image.save("test.pdb", -2);

    // int i = 0;
    // for (const auto& im : image.images()) {
    //     PlotImage plot(im);
    //     plot.plot_atoms(-1);
    //     plot.save("temp/" + std::to_string(++i) + ".png");
    // }

    auto header = image.get_header();
    std::cout << *header << std::endl;

    // image.fit("data/A2M_ma.RSR");

    setting::axes::scattering_intensity_plot_binned_width = 1;

    ScatteringHistogram h(image.get_histogram(-2));

    PlotDistance distance(h);
    PlotIntensity intensity(h);
    distance.save("distance.pdf");
    intensity.save("intensity.pdf");
    return 0;
}