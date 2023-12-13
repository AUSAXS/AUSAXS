#include <em/detail/ImageStackBase.h>
#include <em/detail/header/data/DummyData.h>
#include <em/detail/header/HeaderFactory.h>
#include <em/manager/ProteinManagerFactory.h>
#include <em/ObjectBounds3D.h>
#include <em/Image.h>
#include <data/Molecule.h>
#include <fitter/Fit.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <settings/EMSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/GridSettings.h>
#include <utility/Exceptions.h>
#include <constants/Constants.h>
#include <utility/Axis3D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

#include <fstream>
#include <cstdint>
#include <numeric>
#include <functional>
#include <memory>

using namespace em;

ImageStackBase::ImageStackBase(const std::vector<Image>& images) : size_x(images[0].N), size_y(images[0].M), size_z(static_cast<unsigned int>(images.size())) {    
    data = images;
    for (unsigned int z = 0; z < size_z; ++z) {
        if (image(z).N != size_x || image(z).M != size_y) {throw except::invalid_argument("ImageStackBase::ImageStackBase: All images must have the same dimensions.");}
        image(z).set_z(z);
    }
    phm = factory::create_manager(this);
}

ImageStackBase::ImageStackBase(const io::ExistingFile& file) {
    constants::filetypes::em_map.validate(file);
    header = em::detail::factory::create_header(file);

    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("ImageStackBase::ImageStackBase: Could not open file \"" + file + "\"");}
    input.read(reinterpret_cast<char*>(header->get_data()), header->get_header_size());

    auto map_axes = header->get_axes();
    size_x = map_axes.x.bins;
    size_y = map_axes.y.bins;
    size_z = map_axes.z.bins;

    read(input);
    phm = factory::create_manager(this);
}

ImageStackBase::~ImageStackBase() = default;

void ImageStackBase::save(double cutoff, const io::File& path) const {
    auto protein = phm->get_protein(cutoff);
    protein->save(path);
}

Image& ImageStackBase::image(unsigned int layer) {return data[layer];}

const Image& ImageStackBase::image(unsigned int layer) const {return data[layer];}

unsigned int ImageStackBase::size() const {return size_z;}

const std::vector<Image>& ImageStackBase::images() const {return data;}

std::unique_ptr<hist::ICompositeDistanceHistogram> ImageStackBase::get_histogram(double cutoff) const {
    return phm->get_histogram(cutoff);
}

std::unique_ptr<hist::ICompositeDistanceHistogram> ImageStackBase::get_histogram(const std::shared_ptr<fitter::EMFit> res) const {
    return get_histogram(res->get_parameter("cutoff").value);
}

observer_ptr<data::Molecule> ImageStackBase::get_protein(double cutoff) const {
    return phm->get_protein(cutoff);
}

unsigned int ImageStackBase::count_voxels(double cutoff) const {
    return std::accumulate(data.begin(), data.end(), 0u, [&cutoff] (unsigned int sum, const Image& im) {return sum + im.count_voxels(cutoff);});
}

template<numeric T>
float read_helper(std::ifstream& istream, unsigned int readsize) {
    T value;
    istream.read(reinterpret_cast<char*>(&value), readsize);
    return value;
}

std::function<float(std::ifstream&, unsigned int)> get_read_function(em::detail::header::DataType data_type) {
    switch (data_type) {
        case em::detail::header::DataType::int8: return read_helper<int8_t>;
        case em::detail::header::DataType::int16: return read_helper<int16_t>;
        case em::detail::header::DataType::uint8: return read_helper<uint8_t>;
        case em::detail::header::DataType::uint16: return read_helper<uint16_t>;
        case em::detail::header::DataType::float16: return read_helper<float>;
        case em::detail::header::DataType::float32: return read_helper<float>;
        default: throw except::invalid_argument("ImageStackBase::get_read_function: Invalid data type");
    }
}

void ImageStackBase::read(std::ifstream& istream) {
    data = std::vector<Image>(size_z, Image(header.get()));
    auto[col, row, sec] = header->get_axis_order();

    // the data is stored in the order of column, row, section
    // we have to convert this format to (x, y, z)
    // first determine the limits of each axis
    unsigned int xm, ym, zm;
    auto set_size = [this] (int axis) {
        switch (axis) {
            case 1: return size_x;
            case 2: return size_y;
            case 3: return size_z;
            default: throw except::invalid_argument("ImageStackBase::read: Invalid axis");
        }
    };

    // set the limits of each axis
    xm = set_size(col);
    ym = set_size(row);
    zm = set_size(sec);

    // define an index array to contain the current indices of each axis
    std::array<unsigned int, 3> i = {0, 0, 0};

    // define a permutated reference to each index 
    unsigned int &x = i[col-1];
    unsigned int &y = i[row-1];
    unsigned int &z = i[sec-1];

    // do the actual reading. Note that the default order is 123, so we have to iterate over z first, then y, then x
    auto readfunc = get_read_function(header->get_data_type());
    for (i[2] = 0; i[2] < zm; i[2]++) {
        for (i[1] = 0; i[1] < ym; i[1]++) {
            for (i[0] = 0; i[0] < xm; i[0]++) {
                index(x, y, z) = readfunc(istream, header->get_byte_size());
            }
        }
    }
    // check that we have read the correct number of bytes
    if (istream.peek() != EOF) {throw except::io_error("ImageStackBase::read: File is larger than expected.");}

    // set z values
    for (unsigned int z = 0; z < size_z; z++) {
        image(z).set_z(z);
    }

    // update dummy volumes so they can cover the map interior
    auto axes = header->get_axes();
    double xwidth = axes.x.width();
    double ywidth = axes.y.width();
    double zwidth = axes.z.width();
    double minwidth = std::min({xwidth, ywidth, zwidth})*settings::em::sample_frequency;
    constants::radius::set_dummy_radius(std::sqrt(2*std::pow(minwidth, 2))/2 + settings::grid::width);
}

float& ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) {
    return data[layer].index(x, y);
}

float ImageStackBase::index(unsigned int x, unsigned int y, unsigned int layer) const {
    return data[layer].index(x, y);
}

observer_ptr<em::detail::header::MapHeader> ImageStackBase::get_header() const {
    return header.get();
}

void ImageStackBase::set_header(std::unique_ptr<em::detail::header::MapHeader> header) {
    this->header = std::move(header);
    for (auto& image : data) {
        image.set_header(this->header.get());
    }
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
    if (_rms == 0) {
        double sum = std::accumulate(data.begin(), data.end(), 0.0, [] (double sum, const Image& image) {return sum + image.squared_sum();});
        _rms = std::sqrt(sum/(size_x*size_y*size_z));
    }
    return _rms;
}

observer_ptr<managers::ProteinManager> ImageStackBase::get_protein_manager() const {
    return phm.get();
}