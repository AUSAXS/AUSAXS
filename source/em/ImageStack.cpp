#include <algorithm>
#include <fstream>

//! ##### Remove ##### !
#include <chrono>
using namespace std::chrono;
//! ################## !

#include <em/NoCulling.h>
#include <em/CounterCulling.h>
#include <em/ImageStack.h>
#include <em/PartialHistogramManager.h>
#include <data/Atom.h>
#include <data/Protein.h>
#include <fitter/SimpleIntensityFitter.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotIntensityFitResiduals.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <minimizer/Golden.h>

#include <filesystem>

using namespace setting::em;
using namespace em;

ImageStack::ImageStack(const vector<Image>& images, unsigned int resolution, CullingStrategyChoice csc) 
    : resolution(resolution), size_x(images[0].N), size_y(images[0].M), size_z(images.size()), phm(std::make_unique<em::PartialHistogramManager>(*this)) {
    
    data = images;
    setup(csc);
}

void ImageStack::validate_extension(string file) const {
    string extension = std::filesystem::path(file).extension().string();
    if (extension == ".ccp4") {return;}
    if (extension == ".map") {return;}
    throw except::invalid_extension("Error in Imagestack::validate_extension: Invalid extension \"" + extension + "\".");
}

ImageStack::ImageStack(string file, unsigned int resolution, CullingStrategyChoice csc) 
    : filename(file), header(std::make_shared<ccp4::Header>()), resolution(resolution), phm(std::make_unique<em::PartialHistogramManager>(*this)) {

    validate_extension(file);

    std::ifstream input(file, std::ios::binary);
    if (!input.is_open()) {throw except::io_error("Error in ImageStack::ImageStack: Could not open file \"" + file + "\"");}

    input.read(reinterpret_cast<char*>(header.get()), sizeof(*header));
    read(input, get_byte_size());
    setup(csc);
}

ImageStack::~ImageStack() = default;

void ImageStack::setup(CullingStrategyChoice csc) {
    determine_staining();
    if (negatively_stained()) {std::transform(charge_levels.begin(), charge_levels.end(), charge_levels.begin(), std::negate<double>());}

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
    std::shared_ptr<Protein> protein = phm->get_protein(cutoff);
    protein->save(path);
}

Image& ImageStack::image(unsigned int layer) {return data[layer];}

const Image& ImageStack::image(unsigned int layer) const {return data[layer];}

size_t ImageStack::size() const {return size_z;}

const vector<Image>& ImageStack::images() const {return data;}

std::unique_ptr<Grid> ImageStack::create_grid(double) const {
    throw except::unexpected("Error in Imagestack::create_grid: Not implemented yet.");
}

hist::ScatteringHistogram ImageStack::get_histogram(double cutoff) const {
    return phm->get_histogram(cutoff);
}

hist::ScatteringHistogram ImageStack::get_histogram(const std::shared_ptr<EMFit> res) const {
    return get_histogram(res->get_parameter("cutoff").value);
}

std::shared_ptr<Protein> ImageStack::get_protein(double cutoff) const {
    return phm->get_protein(cutoff);
}

void ImageStack::set_staining(double val) {
    if (val == 0) {throw except::invalid_argument("Error in ImageStack::set_staining: Argument is 0, sign cannot be deduced.");}
    staining = val < 0 ? Staining::negative : Staining::positive;
}

void ImageStack::set_staining(Staining staining) noexcept {
    this->staining = staining;
}

void ImageStack::update_charge_levels(Limit limit) const noexcept {
    setting::em::charge_levels.clear();
    for (unsigned int i = 0; i < 20; i++) {
        setting::em::charge_levels.push_back(limit.min + i*limit.span()/20);
    }
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(const hist::ScatteringHistogram& h) {
    Limit lim = {0, header->dmax};
    mini::Parameter param("cutoff", lim.center(), lim);
    fit(h, param);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(const hist::ScatteringHistogram& h, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(h);} // ensure parameter bounds are present
    SimpleIntensityFitter fitter(h, get_limits());
    return fit_helper(fitter, param);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(const hist::ScatteringHistogram& h) {
    Limit lim = {0, header->dmax};
    mini::Parameter param("cutoff", lim.center(), lim);
    fit(h, param);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit(string filename, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(filename);} // ensure parameter bounds are present
    SimpleIntensityFitter fitter(filename);
    return fit_helper(fitter, param);
}

std::shared_ptr<ImageStack::EMFit> ImageStack::fit_helper(SimpleIntensityFitter& fitter, mini::Parameter param) {
    update_charge_levels(*param.bounds);
    determine_minimum_bounds(param.bounds->min);
    auto func = prepare_function(fitter);

    mini::Golden minimizer(func, param);
    auto res = minimizer.minimize();

    std::shared_ptr<EMFit> emfit = std::make_shared<EMFit>(fitter, res, res.fval);
    emfit->evaluated_points = minimizer.get_evaluated_points();
    return emfit;
}

std::function<double(const double*)> ImageStack::prepare_function(SimpleIntensityFitter& fitter) {
    // convert the calculated intensities to absolute scale
    utility::print_warning("Warning in ImageStack::prepare_function: Not using absolute scale.");
    // auto protein = phm->get_protein(1);
    // double c = setting::em::concentration;                                // concentration
    // double m = protein->get_absolute_mass()*constants::unit::mg;          // mass
    // double DrhoV2 = std::pow(protein->get_relative_charge(), 2);          // charge
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2); // squared scattering length
    // double I0 = DrhoV2*re2*c/m;
    // fitter.normalize_intensity(I0);

    // fit function
    static unsigned int counter;
    counter = 0; // must be in separate line since we want to reset it every time this function is called
    std::function<double(const double*)> chi2 = [&] (const double* params) {
        fitter.set_scattering_hist(get_histogram(params[0]));
        double val = fitter.fit()->fval;
        std::cout << "Step " << counter++ << ": Evaluated cutoff value " << params[0] << " with chi2 " << val << std::endl;
        return val;
    }; 
    return chi2;
}

Dataset ImageStack::cutoff_scan(const Axis& points, string file) {
    SimpleIntensityFitter fitter(file);
    return cutoff_scan_helper(points, fitter);
}

Dataset ImageStack::cutoff_scan(const Axis& points, const hist::ScatteringHistogram& h) {
    SimpleIntensityFitter fitter(h, get_limits());
    return cutoff_scan_helper(points, fitter);
}

Dataset ImageStack::cutoff_scan_helper(const Axis& points, SimpleIntensityFitter& fitter) {
    determine_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    return minimizer.landscape(points.bins);
}

ImageStack::Landscape ImageStack::cutoff_scan_fit(const Axis& points, const hist::ScatteringHistogram& h) {
    SimpleIntensityFitter fitter(h, get_limits());
    determine_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    // cutoff scan
    Landscape landscape;
    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    landscape.contour = minimizer.landscape(points.bins);

    // fit
    double guess; Limit bounds;
    if (positively_stained()) {guess = 2; bounds = Limit(1, 10);}   //! Find a better approach - this will not always work
    else {guess = -2; bounds = Limit(-10, -1);}                     //! Same
    minimizer.clear_parameters();
    minimizer.add_parameter({"cutoff", guess, bounds});
    auto res = minimizer.minimize();

    EMFit emfit(fitter, res, res.fval);
    emfit.evaluated_points = minimizer.get_evaluated_points();
    landscape.fit = emfit;

    return landscape;
}

unsigned int ImageStack::count_voxels(double cutoff) const {
    return std::accumulate(data.begin(), data.end(), 0, [&cutoff] (double sum, const Image& im) {return sum + im.count_voxels(cutoff);});
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

bool ImageStack::positively_stained() const {return staining == Staining::positive;}

bool ImageStack::negatively_stained() const {return staining == Staining::negative;}

void ImageStack::determine_staining() {
    staining = Staining::positive;
    // // we count how many images where the maximum density is positive versus negative
    // double sign = 0;
    // for (unsigned int z = 0; z < size_z; z++) {
    //     Limit limit = image(z).limits();
    //     double min = std::abs(limit.min), max = std::abs(limit.max);
    //     if (1 <= min && max+1 < min) {
    //         sign--;
    //     } else if (1 <= max && min+1 < max) {
    //         sign++;
    //     }
    // }

    // // set the staining type so we don't have to calculate it again later
    // if (sign > 0) {staining = Staining::positive;}
    // else {staining = Staining::negative;}
}

ObjectBounds3D ImageStack::minimum_volume(double cutoff) {
    ObjectBounds3D bounds(size_x, size_y, size_z);
    for (unsigned int z = 0; z < size_z; z++) {
        bounds[z] = image(z).setup_bounds(cutoff);
    }

    return bounds;
}

void ImageStack::determine_minimum_bounds(double min_val) {
    set_staining(min_val);
    double cutoff = positively_stained() ? std::abs(min_val) : -std::abs(min_val);
    std::for_each(data.begin(), data.end(), [&cutoff] (Image& image) {image.setup_bounds(cutoff);});
}