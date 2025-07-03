// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <plots/PlotImage.h>
#include <settings/PlotSettings.h>
#include <em/Image.h>
#include <dataset/SimpleDataset.h>
#include <hist/Histogram2D.h>

using namespace ausaxs::plots;

PlotImage::PlotImage(const em::Image& image, const PlotOptions& options) {
    plot(image, options);
}

PlotImage::~PlotImage() = default;

PlotImage& PlotImage::plot_atoms(const em::Image& image, double cutoff) {
    std::list<data::EMAtom> atoms = image.generate_atoms(cutoff);
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(atoms.size());
    y.reserve(atoms.size());
    for (const data::EMAtom& atom : atoms) {
        x.push_back(atom.coordinates().x());
        y.push_back(atom.coordinates().y());
    }

    SimpleDataset p(x, y);
    ss << "PlotImageAtoms"
        << p.to_string()
        << "\n"
        << plots::PlotOptions("points", {{"color", style::color::black}, {"marker_size", 2}}).to_string()
        << std::endl;
    return *this;
}

void PlotImage::plot(const em::Image& image, const PlotOptions& options) {
    hist::Histogram2D h = image.as_hist();

    auto opt = options;
    opt.draw_markers = true;
    opt.draw_line = false;
    opt.xlabel = "Length [Å]";
    opt.ylabel = "Length [Å]";
    opt.zlabel = "Electron Density [Arb.]";

    std::string contours;
    for (auto c : settings::plots::contour) {
        contours += std::to_string(c) + " ";
    }

    ss << "PlotImage\n"
        << h.to_string()
        << "\n"
        << options.to_string()
        << contours
        << std::endl;
}

void PlotImage::quick_plot(const em::Image& image, const PlotOptions& options, const io::File& path) {
    PlotImage p(image, options);
    p.save(path);
}