#include <memory>
#include <list>
#include <algorithm>
#include <filesystem>
#include <cmath>

//! ##### Remove ##### !
#include <chrono>
using namespace std::chrono;
//! ################## !

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

ImageStack::ImageStack(const vector<Image>& images, unsigned int resolution, CullingStrategyChoice csc) : size_x(images[0].N), size_y(images[0].M), size_z(images.size()) {
    data = images;
    setup(csc);
}

ImageStack::ImageStack(string file, unsigned int resolution, CullingStrategyChoice csc) : filename(file), header(std::make_shared<ccp4::Header>()), resolution(resolution) {
    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("Error in ImageStack::ImageStack: Could not open file \"" + file + "\"");}

    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    read(input, get_byte_size());
    setup(csc);
}

void ImageStack::setup(CullingStrategyChoice csc) {
    determine_staining();

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

size_t ImageStack::size() const {return size_z;}

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

std::unique_ptr<Grid> ImageStack::create_grid(double) const {
    std::cout << "Error in ImageStack::create_grid: Not implemented yet. " << std::endl;
    exit(1);
}

ScatteringHistogram ImageStack::get_histogram(double cutoff) const {
    // auto start = high_resolution_clock::now();
    std::unique_ptr<Protein> protein = create_protein(cutoff);
    // auto duration = duration_cast<seconds>(high_resolution_clock::now() - start);
    // std::cout << "\tTook " << duration.count() << " seconds to prepare protein with " << protein->atom_size() << " atoms." << std::endl;
    // protein->generate_new_hydration();
    return protein->get_histogram();
}

ScatteringHistogram ImageStack::get_histogram(const std::shared_ptr<EMFit> res) const {
    return get_histogram(res->params.at("cutoff"));
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(const ScatteringHistogram& h) {
    SimpleIntensityFitter fitter(h, get_limits());
    determine_minimum_bounds();
    return fit_helper(fitter);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(string filename) {
    SimpleIntensityFitter fitter(filename);
    determine_minimum_bounds();
    return fit_helper(fitter);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit_helper(SimpleIntensityFitter& fitter) {
    // fit function
    unsigned int counter = 0;
    std::function<double(const double*)> chi2 = [&] (const double* params) {
        std::cout << "\nStep " << ++counter << std::endl;
        auto start = high_resolution_clock::now();
        fitter.set_scattering_hist(get_histogram(params[0]));
        double val = fitter.fit()->chi2;
        auto duration = duration_cast<seconds>(high_resolution_clock::now() - start);
        std::cout << "\tTotal time: " << duration.count() << std::endl;

        std::cout << "\tEvaluated cutoff value " << params[0] << " with chi2 " << val << std::endl;
        return val;
    }; 

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
    chi2(result);

    return std::make_shared<EMFit>(fitter, minimizer, minimizer->MinValue());
}

size_t ImageStack::get_byte_size() const {
    return header->get_byte_size();
}

void ImageStack::read(std::ifstream& istream, size_t byte_size) {
    size_x = header->nx; size_y = header->ny; size_z = header->nz;

    data = vector<Image>(size_z, Image(header));
    for (unsigned int i = 0; i < size_x; i++) {
        for (unsigned int j = 0; j < size_y; j++) {
            for (unsigned int k = 0; k < size_z; k++) {
                istream.read(reinterpret_cast<char*>(&index(i, j, k)), byte_size);
            }
        }
    }

    // set z values
    for (unsigned int z = 0; z < size_z; z++) {
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
    return resolution == 0 ? Limit(setting::fit::q_low, setting::fit::q_high) : Limit(setting::fit::q_low, 2*M_PI/resolution);
}

double ImageStack::mean() const {
    double sum = 0;
    for (unsigned int z = 0; z < size_z; z++) {
        sum += image(z).mean();
    }
    return sum/size_z;
}

bool ImageStack::is_positively_stained() const {return staining > 0;}

void ImageStack::determine_staining() {
    // we count how many images where the maximum density is positive versus negative
    double sign = 0;
    for (unsigned int z = 0; z < size_z; z++) {
        Limit limit = image(z).limits();
        double min = std::abs(limit.min), max = std::abs(limit.max);
        if (1 <= min && max+1 < min) {
            sign--;
        } else if (1 <= max && min+1 < max) {
            sign++;
        }
    }

    // set the staining type so we don't have to calculate it again later
    if (sign > 0) {staining = 1;}
    else {staining = -1;}
}

ObjectBounds3D ImageStack::minimum_volume(double cutoff) {
    ObjectBounds3D bounds(size_x, size_y, size_z);
    for (unsigned int z = 0; z < size_z; z++) {
        bounds[z] = image(z).setup_bounds(cutoff);
    }

    return bounds;
}

void ImageStack::determine_minimum_bounds() {
    double cutoff = is_positively_stained() ? 1 : -1;
    std::for_each(data.begin(), data.end(), [&cutoff] (Image& image) {image.setup_bounds(cutoff);});
}