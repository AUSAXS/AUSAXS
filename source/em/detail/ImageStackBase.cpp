#include <algorithm>
#include <fstream>
#include <filesystem>

#include <preprocessor.h>
#include <em/NoCulling.h>
#include <em/CounterCulling.h>
#include <em/detail/ImageStackBase.h>
#include <em/PartialHistogramManager.h>
#include <data/Atom.h>
#include <data/Protein.h>
#include <fitter/SimpleIntensityFitter.h>
#include <plots/all.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <mini/all.h>
#include <math/Statistics.h>

using namespace em;

ImageStackBase::ImageStackBase(const std::vector<Image>& images) : size_x(images[0].N), size_y(images[0].M), size_z(images.size()), phm(std::make_unique<em::PartialHistogramManager>(*this)) {    
    data = images;
    phm->set_charge_levels(); // set default charge levels
}

ImageStackBase::ImageStackBase(std::string file) : filename(file), header(std::make_shared<ccp4::Header>()), phm(std::make_unique<em::PartialHistogramManager>(*this)) {
    constants::filetypes::em_map.validate(file);

    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("ImageStackBase::ImageStackBase: Could not open file \"" + file + "\"");}

    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    size_x = header->nx; size_y = header->ny; size_z = header->nz;
    read(input, get_byte_size());
    phm->set_charge_levels(); // set default charge levels
}

ImageStackBase::~ImageStackBase() = default;

void ImageStackBase::save(double cutoff, std::string path) const {
    std::shared_ptr<Protein> protein = phm->get_protein(cutoff);
    protein->save(path);
}

Image& ImageStackBase::image(unsigned int layer) {return data[layer];}

const Image& ImageStackBase::image(unsigned int layer) const {return data[layer];}

unsigned int ImageStackBase::size() const {return size_z;}

const std::vector<Image>& ImageStackBase::images() const {return data;}

hist::ScatteringHistogram ImageStackBase::get_histogram(double cutoff) const {
    return phm->get_histogram(cutoff);
}

hist::ScatteringHistogram ImageStackBase::get_histogram(const std::shared_ptr<EMFit> res) const {
    return get_histogram(res->get_parameter("cutoff").value);
}

std::shared_ptr<Protein> ImageStackBase::get_protein(double cutoff) const {
    return phm->get_protein(cutoff);
}

unsigned int ImageStackBase::count_voxels(double cutoff) const {
    return std::accumulate(data.begin(), data.end(), 0, [&cutoff] (double sum, const Image& im) {return sum + im.count_voxels(cutoff);});
}

size_t ImageStackBase::get_byte_size() const {
    return header->get_byte_size();
}

void ImageStackBase::read(std::ifstream& istream, size_t byte_size) {
    data = std::vector<Image>(size_z, Image(header));

    int col = header->mapc; // column axis
    int row = header->mapr; // row axis
    int sec = header->maps; // section axis

    unsigned int i1, i2, i3;
    auto set_size = [this] (int axis) {
        switch (axis) {
            case 1: return size_x;
            case 2: return size_y;
            case 3: return size_z;
            default: throw except::invalid_argument("ImageStackBase::read: Invalid axis");
        }
    };

    i1 = set_size(col);
    i2 = set_size(row);
    i3 = set_size(sec);

    for (unsigned int i = 0; i < i3; i++) {
        for (unsigned int j = 0; j < i2; j++) {
            for (unsigned int k = 0; k < i1; k++) {
                istream.read(reinterpret_cast<char*>(&index(i, j, k)), byte_size);
            }
        }
    }

    // for (unsigned int i = 0; i < size_x; i++) {
    //     for (unsigned int j = 0; j < size_y; j++) {
    //         for (unsigned int k = 0; k < size_z; k++) {
    //             istream.read(reinterpret_cast<char*>(&index(i, j, k)), byte_size);
    //         }
    //     }
    // }
    if (istream.peek() != EOF) {throw except::io_error("ImageStackBase::read: File is larger than expected.");}

    // set z values
    for (unsigned int z = 0; z < size_z; z++) {
        image(z).set_z(z);
    }
}

float& ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) {
    return data[layer].index(x, y);
}

float ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) const {
    return data[layer].index(x, y);
}

std::shared_ptr<ccp4::Header> ImageStackBase::get_header() const {
    return header;
}

Limit ImageStackBase::get_limits() const {
    return Limit(setting::axes::qmin, setting::axes::qmax);
}

double ImageStackBase::mean() const {
    double sum = 0;
    for (unsigned int z = 0; z < size_z; z++) {
        sum += image(z).mean();
    }
    return sum/size_z;
}

ObjectBounds3D ImageStackBase::minimum_volume(double cutoff) {
    ObjectBounds3D bounds(size_x, size_y, size_z);
    for (unsigned int z = 0; z < size_z; z++) {
        bounds[z] = image(z).setup_bounds(cutoff);
    }

    return bounds;
}

void ImageStackBase::set_minimum_bounds(double min_val) {
    std::for_each(data.begin(), data.end(), [&min_val] (Image& image) {image.setup_bounds(min_val);});
}

double ImageStackBase::from_level(double sigma) const {
    return sigma*rms();
}

double ImageStackBase::to_level(double cutoff) const {
    return cutoff/rms();
}

double ImageStackBase::rms() const {
    static double rms = 0; // only initialized once
    if (rms == 0) {
        double sum = std::accumulate(data.begin(), data.end(), 0.0, [] (double sum, const Image& image) {return sum + image.squared_sum();});
        rms = std::sqrt(sum/(size_x*size_y*size_z));
    }
    return rms;
}

std::shared_ptr<em::PartialHistogramManager> ImageStackBase::get_histogram_manager() const {
    return phm;
}