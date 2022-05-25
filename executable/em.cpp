#include <iostream>

#include <em/ImageStack.h>
#include <plots/PlotImage.h>
#include <plots/PlotIntensity.h>
#include <plots/PlotDistance.h>
#include <utility/Utility.h>

using std::string;

int main(int argc, char const *argv[]) {
    // em::ImageStack image("data/A2M_map.ccp4");
    // image.plot(std::stoi(argv[1]));
    // image.fit("data/A2M_ma.RSR");

    // string file = "data/A2M_ma.ccp4";
    string file = "sim/native_25.ccp4";
    em::ImageStack image(file); 
    std::cout << image.get_header()->to_string() << std::endl;
    // plots::PlotImage plot(image.image(std::stoi(argv[1])));
    // plot.plot_atoms(-1);
    // plot.save("temp.pdf");
    // image.save("test.pdb", -2);

    exit(1);

    int i = 0;
    for (const auto& im : image.images()) {
        plots::PlotImage plot(im);
        // plot.plot_atoms(-1);
        plot.save("figures/em/" + utility::stem(file) + "/" + std::to_string(++i) + ".png");
    }

    // setting::axes::scattering_intensity_plot_binned_width = 0.01;
    // setting::protein::use_effective_charge = false;

    // image.fit("data/A2M_ma.RSR");

    // ScatteringHistogram h(image.get_histogram(std::atof(argv[1])));
    // plots::PlotDistance distance(h);
    // plots::PlotIntensity intensity(h);
    // distance.save("distance.pdf");
    // intensity.save("intensity.pdf");
    return 0;
}