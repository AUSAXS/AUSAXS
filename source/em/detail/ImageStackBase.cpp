#include <em/detail/ImageStackBase.h>
#include <em/manager/ProteinManagerFactory.h>
#include <em/ObjectBounds3D.h>
#include <em/Image.h>
#include <data/Protein.h>
#include <fitter/Fit.h>
#include <mini/detail/FittedParameter.h>
#include <settings/EMSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/Exceptions.h>
#include <utility/Constants.h>
#include <mini/detail/Evaluation.h>
#include <Symbols.h>

#include <fstream>

em::ImageStackBase::ImageStackBase(const std::vector<Image>& images) : size_x(images[0].N), size_y(images[0].M), size_z(images.size()) {    
    data = images;
    phm = em::factory::create_manager(this);
}

em::ImageStackBase::ImageStackBase(const io::ExistingFile& file) : filename(file), header(std::make_shared<ccp4::Header>()) {
    constants::filetypes::em_map.validate(file);

    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("ImageStackBase::ImageStackBase: Could not open file \"" + file + "\"");}
    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    size_x = header->nx; size_y = header->ny; size_z = header->nz;
    read(input, get_byte_size());
    phm = em::factory::create_manager(this);
}

em::ImageStackBase::~ImageStackBase() = default;

void em::ImageStackBase::save(double cutoff, std::string path) const {
    std::shared_ptr<Protein> protein = phm->get_protein(cutoff);
    protein->save(path);
}

em::Image& em::ImageStackBase::image(unsigned int layer) {return data[layer];}

const em::Image& em::ImageStackBase::image(unsigned int layer) const {return data[layer];}

unsigned int em::ImageStackBase::size() const {return size_z;}

const std::vector<em::Image>& em::ImageStackBase::images() const {return data;}

hist::ScatteringHistogram em::ImageStackBase::get_histogram(double cutoff) const {
    return phm->get_histogram(cutoff);
}

hist::ScatteringHistogram em::ImageStackBase::get_histogram(const std::shared_ptr<fitter::EMFit> res) const {
    return get_histogram(res->get_parameter("cutoff").value);
}

std::shared_ptr<Protein> em::ImageStackBase::get_protein(double cutoff) const {
    std::cout << "phm is " << (phm == nullptr ? "null" : "not null") << std::endl;
    return phm->get_protein(cutoff);
}

unsigned int em::ImageStackBase::count_voxels(double cutoff) const {
    return std::accumulate(data.begin(), data.end(), 0, [&cutoff] (double sum, const Image& im) {return sum + im.count_voxels(cutoff);});
}

unsigned int em::ImageStackBase::get_byte_size() const {
    return header->get_byte_size();
}

void em::ImageStackBase::read(std::ifstream& istream, unsigned int byte_size) {
    data = std::vector<Image>(size_z, Image(header));

    int col = header->mapc; // column axis
    int row = header->mapr; // row axis
    int sec = header->maps; // section axis

    // the data is stored in the order of column, row, section
    // we have to convert this format to (x, y, z)
    // first determine the limits of each axis
    unsigned int xm, ym, zm;
    auto set_size = [this] (int axis) {
        switch (axis) {
            case 1: return size_x;
            case 2: return size_y;
            case 3: return size_z;
            default: throw except::invalid_argument("ImageStack::read: Invalid axis");
        }
    };

    // set the limits of each axis
    xm = set_size(col);
    ym = set_size(row);
    zm = set_size(sec);

    // define an index array to contain the current indices of each axis
    std::array<unsigned int, 3> i = {0, 0, 0};

    // define a permutated reference to each index 
    unsigned int &x = i[header->mapc-1];
    unsigned int &y = i[header->mapr-1];
    unsigned int &z = i[header->maps-1];

    // do the actual reading. Note that the default order is 123, so we have to iterate over z first, then y, then x
    for (i[2] = 0; i[2] < zm; i[2]++) {
        for (i[1] = 0; i[1] < ym; i[1]++) {
            for (i[0] = 0; i[0] < xm; i[0]++) {
                istream.read(reinterpret_cast<char*>(&index(x, y, z)), byte_size);
            }
        }
    }
    // check that we have read the correct number of bytes
    if (istream.peek() != EOF) {throw except::io_error("ImageStack::read: File is larger than expected.");}

    // set z values
    for (unsigned int z = 0; z < size_z; z++) {
        image(z).set_z(z);
    }
}

float& em::ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) {
    return data[layer].index(x, y);
}

float em::ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) const {
    return data[layer].index(x, y);
}

std::shared_ptr<em::ccp4::Header> em::ImageStackBase::get_header() const {
    return header;
}

void em::ImageStackBase::set_header(std::shared_ptr<ccp4::Header> header) {
    this->header = header;
    for (auto& image : data) {
        image.set_header(header);
    }
}

Limit em::ImageStackBase::get_limits() const {
    return Limit(settings::axes::qmin, settings::axes::qmax);
}

double em::ImageStackBase::mean() const {
    double sum = 0;
    for (unsigned int z = 0; z < size_z; z++) {
        sum += image(z).mean();
    }
    return sum/size_z;
}

em::ObjectBounds3D em::ImageStackBase::minimum_volume(double cutoff) {
    ObjectBounds3D bounds(size_x, size_y, size_z);
    for (unsigned int z = 0; z < size_z; z++) {
        bounds[z] = image(z).setup_bounds(cutoff);
    }

    return bounds;
}

void em::ImageStackBase::set_minimum_bounds(double min_val) {
    std::for_each(data.begin(), data.end(), [&min_val] (Image& image) {image.setup_bounds(min_val);});
}

double em::ImageStackBase::from_level(double sigma) const {
    return sigma*rms();
}

double em::ImageStackBase::to_level(double cutoff) const {
    return cutoff/rms();
}

double em::ImageStackBase::rms() const {
    static double rms = 0; // only initialized once
    if (rms == 0) {
        double sum = std::accumulate(data.begin(), data.end(), 0.0, [] (double sum, const Image& image) {return sum + image.squared_sum();});
        rms = std::sqrt(sum/(size_x*size_y*size_z));
    }
    return rms;
}

std::shared_ptr<em::managers::ProteinManager> em::ImageStackBase::get_protein_manager() const {
    return phm;
}