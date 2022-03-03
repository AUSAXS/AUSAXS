#include <memory>
#include <list>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <em/image.h>
#include <data/Atom.h>
#include <data/Protein.h>
#include <fitter/SimpleIntensityFitter.h>
#include "plots/PlotIntensityFit.h"
#include "plots/PlotIntensityFitResiduals.h"
#include <Exceptions.h>

using namespace em;
using std::list;

ImageStack::ImageStack(string file) : header(std::make_shared<ccp4::Header>()) {
    std::ifstream input("data/A2M_map.ccp4", std::ios::binary);
    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    read(input, get_byte_size());
}

Image& ImageStack::image(unsigned int layer) {return data[layer];}

size_t ImageStack::size() const {return header->nz;}

const vector<Image>& ImageStack::images() const {return data;}

ImageStack::ImageStack(std::shared_ptr<ccp4::Header> header, std::ifstream& istream) : header(header) {
    read(istream, get_byte_size());
}

std::unique_ptr<Protein> ImageStack::create_protein(double cutoff) const {
    // we use a list since we will have to append quite a few other lists to it
    list<Atom> atoms;
    int i = 0;
    for (const Image& image: data) {
        list<Atom> im_atoms = image.generate_atoms(cutoff);
        std::cout << "Generating " << im_atoms.size() << " atoms for image " << ++i << std::endl;
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
    }

    // convert list to vector
    vector<Atom> v;
    v.reserve(atoms.size());
    v.assign(std::move_iterator(atoms.begin()), std::move_iterator(atoms.end()));
    return std::make_unique<Protein>(v);
}

std::unique_ptr<Grid> ImageStack::create_grid(double cutoff) const {
    vector<Atom> atoms;
    atoms.reserve(header->nx);

    double xscale = header->cella_x/header->nx;
    double yscale = header->cella_y/header->ny;
    double zscale = header->cella_z/header->nz;
    for (int x = 0; x < header->nx; x++) {
        for (int y = 0; y < header->ny; y++) {
            for (int z = 0; z < header->nz; z++) {
                float val = index(x, y, z);
                if (val < cutoff) {
                    continue;
                }

                Vector3 coords{x*xscale, y*yscale, z*zscale};
                atoms.push_back(Atom(coords, val, "C", "C", 0));
            }
        }
    }

    Protein protein(atoms);
    return std::make_unique<Grid>(atoms);
}

ScatteringHistogram ImageStack::calc_scattering_hist() const {
    std::unique_ptr<Protein> protein = create_protein(-1);
    protein->generate_new_hydration();
    return protein->get_histogram();
}

void ImageStack::fit(string filename) const {
    SimpleIntensityFitter fitter(filename, calc_scattering_hist());
    std::shared_ptr<Fitter::Fit> result = fitter.fit();

    // Fit plot
    PlotIntensityFit plot_f(fitter);
    plot_f.save("em_intensity_fit." + setting::figures::format);

    // Residual plot
    PlotIntensityFitResiduals plot_r(fitter);
    plot_r.save("em_residuals." + setting::figures::format);
}

size_t ImageStack::get_byte_size() const {
    return header->get_byte_size();
}

void ImageStack::plot_without_solution(unsigned int layer) const {
    data[layer].plot_without_solution();
}

void ImageStack::plot() const {
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "canvas", 600, 600);
    std::unique_ptr<TH3D> hist = std::make_unique<TH3D>("hist", "hist", header->nx, 0, 1, header->ny, 0, 1, header->nz, 0, 1);
    for (int z = 0; z < header->nz; z++) {
        for (int y = 0; y < header->ny; y++) {
            for (int x = 0; x < header->nx; x++) {
                hist->SetBinContent(x, y, z, index(x, y, z));
            }
        }
    }

    hist->Draw("box");
    canvas->SaveAs("test3d.pdf");
}

void ImageStack::read(std::ifstream& istream, size_t byte_size) {
    data = vector<Image>(header->nz, Image(header));
    for (int i = 0; i < header->nx; i++) {
        for (int j = 0; j < header->ny; j++) {
            for (int k = 0; k < header->nz; k++) {
                istream.read(reinterpret_cast<char*>(&index(i, j, k)), byte_size);
            }
        }
    }
}

float& ImageStack::index(unsigned int x, unsigned int y, unsigned int layer) {
    return data[layer].index(x, y);
}

float ImageStack::index(unsigned int x, unsigned int y, unsigned int layer) const {
    return data[layer].index(x, y);
}

std::shared_ptr<ccp4::Header> ImageStack::get_header() const {
    return header;
}