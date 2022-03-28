#include <memory>
#include <list>
#include <algorithm>
#include <filesystem>

#include <em/NoCulling.h>
#include <em/CounterCulling.h>
#include <em/ImageStack.h>
#include <data/Atom.h>
#include <data/Protein.h>
#include <fitter/SimpleIntensityFitter.h>
#include "plots/PlotIntensityFit.h"
#include "plots/PlotIntensityFitResiduals.h"
#include <Exceptions.h>

using namespace setting::em;
using namespace em;

ImageStack::ImageStack(string file, CullingStrategyChoice csc) : filename(file), header(std::make_shared<ccp4::Header>()) {
    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("Error in ImageStack::ImageStack: Could not open file \"" + file + "\"");}

    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    read(input, get_byte_size());
    setup(csc);
}

void ImageStack::setup(CullingStrategyChoice csc) {
    switch (csc) {
        case CullingStrategyChoice::NoStrategy:
            culler = std::make_unique<NoCulling>();
            break;
        case CullingStrategyChoice::CounterStrategy:
            culler = std::make_unique<CounterCulling>();
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown PlacementStrategy");
    }
}

void ImageStack::save(string path, double cutoff) const {
    std::unique_ptr<Protein> protein = create_protein(cutoff);
    protein->save(path);
}

Image& ImageStack::image(unsigned int layer) {return data[layer];}

const Image& ImageStack::image(unsigned int layer) const {return data[layer];}

size_t ImageStack::size() const {return header->nz;}

const vector<Image>& ImageStack::images() const {return data;}

std::unique_ptr<Protein> ImageStack::create_protein(double cutoff) const {
    // we use a list since we will have to append quite a few other lists to it
    std::list<Atom> atoms;
    for (const Image& image: data) {
        std::list<Atom> im_atoms = image.generate_atoms(cutoff);
        atoms.splice(atoms.end(), im_atoms); // move im_atoms to end of atoms
    }

    // convert list to vector
    vector<Atom> v = culler->cull(atoms);    
    return std::make_unique<Protein>(v);
}

std::unique_ptr<Grid> ImageStack::create_grid(double cutoff) const {
    std::cout << "Error in ImageStack::create_grid: Not implemented yet. " << std::endl;
    exit(1);

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

ScatteringHistogram ImageStack::get_histogram(double cutoff) const {
    std::unique_ptr<Protein> protein = create_protein(cutoff);
    // protein->generate_new_hydration();
    return protein->get_histogram();
}

void ImageStack::fit(const ScatteringHistogram& h) const {
    SimpleIntensityFitter fitter(h, get_limits());
    fit_helper(fitter);
}

void ImageStack::fit(string filename) const {
    SimpleIntensityFitter fitter(filename);
    fit_helper(fitter);
}

void ImageStack::fit_helper(SimpleIntensityFitter& fitter) const {
    // fit function
    unsigned int counter = 0;
    std::function<double(const double*)> chi2 = [&] (const double* params) {
        fitter.set_scattering_hist(get_histogram(params[0]));
        double val = fitter.fit()->chi2;
        std::cout << ++counter << ": " << params[0] << " with chi2: " << val << std::endl;
        return val;
    }; 

    std::cout << "\tMINIMIZER STARTED " << std::endl;
    // perform the fit
    ROOT::Math::Functor functor = ROOT::Math::Functor(chi2, 1);
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "migrad"); 
    minimizer->SetFunction(functor);

    if (is_positively_stained()) {minimizer->SetLimitedVariable(0, "cutoff", 4, 1, 1, 10);}
    else {minimizer->SetLimitedVariable(0, "cutoff", -4, 1, -1, -10);}

    minimizer->SetStrategy(2);
    minimizer->SetPrintLevel(2);
    minimizer->Minimize();
    const double* result = minimizer->X();

    // set optimal cutoff
    std::cout << "Optimal cutoff is " << result[0] << std::endl;
    fitter.set_scattering_hist(get_histogram(result[0]));
    fitter.fit();

    // Fit plot
    plots::PlotIntensityFit plot_f(fitter);
    plot_f.save("em_intensity_fit." + setting::figures::format);

    // Residual plot
    plots::PlotIntensityFitResiduals plot_r(fitter);
    plot_r.save("em_residuals." + setting::figures::format);
}

size_t ImageStack::get_byte_size() const {
    return header->get_byte_size();
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

    // set z values
    for (int z = 0; z < header->nz; z++) {
        image(z).set_z(z);
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

Limit ImageStack::get_limits() const {
    // get stem of filename
    auto p = std::filesystem::path(filename);
    string stem = p.stem();

    // extract any number that is part of the stem
    string num;
    std::copy_if(stem.begin(), stem.end(), std::back_inserter(num), [] (char c) {return std::isdigit(c);});

    // if no number is present, use default q-values
    double max;
    if (num.empty()) {
        std::cout << "Could not parse resolution in filename \"" + filename + "\"; using default limits." << std::endl;
        max = setting::fit::q_high;
    }

    // otherwise we assume it was the resolution
    else {
        std::cout << "Assuming \"" + num + "\" in filename \"" + filename + "\" refers to the resolution." << std::endl;
        max = 2*M_PI/std::stod(num);
    }

    return Limit(setting::fit::q_low, max);
}

double ImageStack::mean() const {
    double sum = 0;
    for (int z = 0; z < header->nz; z++) {
        sum += image(z).mean();
    }
    return sum/header->nz;
}

bool ImageStack::is_positively_stained() const {
    // we count how many images where the maximum density is positive versus negative
    double sign = 0;
    for (int z = 0; z < header->nz; z++) {
        Limit limit = image(z).limits();
        if (std::abs(limit.min) < std::abs(limit.max)) {
            sign--;
        } else {
            sign++;
        }
    }
    return sign > 0;
}