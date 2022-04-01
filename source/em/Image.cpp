#include <list>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <em/Image.h>

using namespace em;
using std::list, std::vector;

Image::Image(std::shared_ptr<ccp4::Header> header, unsigned int layer) : N(header->nx), M(header->ny), header(header), data(N, M), z(layer), bounds(N, M) {}

Image::Image(const Matrix<float>& data) : N(data.N), M(data.M), data(data), bounds(N, M) {}

void Image::set_z(unsigned int z) {
    this->z = z;
}

float Image::index(unsigned int x, unsigned int y) const {return data.index(x, y);}
float& Image::index(unsigned int x, unsigned int y) {return data.index(x, y);}

list<Atom> Image::generate_atoms(double cutoff) const {
    if (header == nullptr) {throw except::invalid_operation("Error in Image::as_hist: Header must be initialized to use this method.");}
    list<Atom> atoms;

    // determine which comparison function to use
    std::function<bool(double)> compare_positive = [&cutoff] (double val) {return val < cutoff;};
    std::function<bool(double)> compare_negative = [&cutoff] (double val) {return val > cutoff;};
    std::function<bool(double)> compare_func = 0 <= cutoff ? compare_positive : compare_negative;

    // loop through all pixels in this image
    double xscale = header->cella_x/N;
    double yscale = header->cella_y/M;
    double zscale = header->cella_z/header->nz;
    for (unsigned int x = 0; x < N; x++) {
        for (unsigned int y = bounds[x].min; y < bounds[x].max; y++) {
            float val = index(x, y);
            if (compare_func(val)) {
                continue;
            }

            Vector3 coords{x*xscale, y*yscale, z*zscale};
            atoms.push_back(Atom(0, "C", "", "LYS", "", 0, "", coords, val, 0, "C", ""));
        }
    }

    return atoms;
}

std::unique_ptr<TH2D> Image::as_hist() const {
    if (header == nullptr) {throw except::invalid_operation("Error in Image::as_hist: Header must be initialized to use this method.");}
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>("hist", "hist", header->nx, 0, header->cella_x, header->ny, 0, header->cella_y);
    for (unsigned int x = 0; x < N; x++) {
        for (unsigned int y = 0; y < M; y++) {
            hist->SetBinContent(x, y, index(x, y));
        }
    }
    return hist;
}

double Image::mean() const {
    double sum = 0;
    for (unsigned int x = 0; x < N; x++) {
        for (unsigned int y = 0; y < M; y++) {
            sum += index(x, y);
        }
    }

    return sum/(N*M);
}

Limit Image::limits() const {
    double min = 1e9, max = -1e9;
    for (unsigned int x = 0; x < N; x++) {
        for (unsigned int y = 0; y < M; y++) {
            double val = index(x, y);
            min = std::min(min, val);
            max = std::max(max, val);
        }
    }

    return Limit(min, max);
}

const ObjectBounds2D& Image::get_bounds() const {
    return bounds;
}

const ObjectBounds2D& Image::setup_bounds(double cutoff) {
    std::function<bool(double)> accept_positive = [&cutoff] (double val) {return val > cutoff;};
    std::function<bool(double)> accept_negative = [&cutoff] (double val) {return val < cutoff;};
    std::function<bool(double)> accept_func = 0 <= cutoff ? accept_positive : accept_negative;

    for (unsigned int x = 0; x < N; x++) {
        bounds[x].min = 0;
        bounds[x].max = 0;
        bool min_set = false;
        for (unsigned int y = 0; y < M; y++) {
            if (!accept_func(index(x, y))) {continue;}
            if (!min_set) {
                bounds[x].min = y; // update min val to this index
                bounds[x].max = y; // also set max in case this is the only entry
                min_set = true;
            } else {
                bounds[x].max = y;
            }
        }
    }

    return bounds;
}