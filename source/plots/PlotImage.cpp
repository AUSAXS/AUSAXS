#include <plots/PlotImage.h>
#include <settings/PlotSettings.h>

plots::PlotImage::PlotImage(const em::Image& image) {
    plot(image);
}

plots::PlotImage::~PlotImage() = default;

void plots::PlotImage::plot_atoms(const em::Image& image, double cutoff) {
    const std::list<Atom>& atoms = image.generate_atoms(cutoff);
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(atoms.size());
    y.reserve(atoms.size());
    for (const Atom& atom : atoms) {
        x.push_back(atom.coords.x());
        y.push_back(atom.coords.y());
    }

    SimpleDataset p(x, y);
    p.add_plot_options("points", {{"color", style::color::black}, {"marker_size", 2}});

    ss << "PlotImageAtoms"
        << p.to_string()
        << "\n"
        << p.get_plot_options().to_string()
        << std::endl;
}

void plots::PlotImage::plot(const em::Image& image) {
    hist::Histogram2D h = image.as_hist();
    h.add_plot_options("points", {{"xlabel", "Length [Å]"}, {"ylabel", "Length [Å]"}, {"zlabel", "Electron Density [Arb.]"}});

    std::string contours;
    for (auto c : settings::plots::contour) {
        contours += std::to_string(c) + " ";
    }

    ss << "PlotImage\n"
        << h.to_string()
        << "\n"
        << h.get_plot_options().to_string()
        << contours
        << std::endl;
}

void plots::PlotImage::quick_plot(const em::Image& image, const io::File& path) {
    plots::PlotImage p(image);
    p.save(path);
}