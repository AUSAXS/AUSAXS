#include <iostream>

#include <em/ImageStack.h>
#include <plots/PlotImage.h>
#include <plots/PlotIntensity.h>
#include <plots/PlotDistance.h>
#include <Exceptions.h>

using std::string;

int main(int argc, char const *argv[]) {
    string pdb_file = argv[1];

    Protein protein(pdb_file);
    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve

    setting::em::max_atoms = 2000;
    gStyle->SetPalette(kSolar);
    auto cols = TColor::GetPalette();

    unsigned int color_step = (cols.GetSize()-1)/argc;
    for (int i = 2; i < argc; i++) {
        std::cout << "Now fitting " << argv[i] << "..." << std::endl;
        em::ImageStack image(argv[i]);
        auto fit = image.fit(pdb_file);
        plot.plot_intensity(image.get_histogram(fit).calc_debye_scattering_intensity(), cols.At(color_step*(i-1)));
    }

    plot.save("figures/stuff.pdf");
    return 0;
}