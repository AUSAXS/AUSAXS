#include <list>
#include <vector>

#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <em/Image.h>
#include <utility/Utility.h>

using namespace em;
using std::list, std::vector;

Image::Image(std::shared_ptr<ccp4::Header> header, unsigned int layer) : N(header->nx), M(header->ny), header(header), data(N, M), z(layer), bounds(N, M) {}

Image::Image(const Matrix<float>& data) : N(data.N), M(data.M), data(data), bounds(N, M) {}

Image::Image(const Matrix<float>& data, std::shared_ptr<ccp4::Header> header, unsigned int layer) : N(data.N), M(data.M), header(header), data(data), z(layer), bounds(N, M) {}

void Image::set_z(unsigned int z) {
    this->z = z;
}

float Image::index(unsigned int x, unsigned int y) const {return data.index(x, y);}
float& Image::index(unsigned int x, unsigned int y) {return data.index(x, y);}

list<Atom> Image::generate_atoms(double cutoff) const {
    if (__builtin_expect(header == nullptr, false)) {throw except::invalid_operation("Error in Image::generate_atoms: Header must be initialized to use this method.");}
    list<Atom> atoms;

    // loop through all pixels in this image
    double xscale = header->cella_x/N;
    double yscale = header->cella_y/M;
    double zscale = header->cella_z/header->nz;

    unsigned int step = setting::em::sample_frequency;
    for (unsigned int x = 0; x < N; x += step) {
        for (unsigned int y = bounds[x].min; y < bounds[x].max; y += step) {
            float val = index(x, y);
            if (val < cutoff) {
                continue;
            }

            Vector3 coords{x*xscale, y*yscale, z*zscale};
            atoms.push_back(Atom(0, "C", "", "LYS", "", 0, "", coords, val, 0, "C", ""));
            // atoms.push_back(Atom(0, "C", "", "LYS", "", 0, "", coords, 1, 0, "C", ""));
        }
    }

    return atoms;
}

unsigned int Image::count_voxels(double cutoff) const {
    list<Atom> l = generate_atoms(cutoff);
    return l.size();
}

double Image::squared_sum() const {
    double sum = 0;
    for (unsigned int x = 0; x < N; x++) {
        for (unsigned int y = 0; y < M; y++) {
            sum += std::pow(index(x, y), 2);
        }
    }
    return sum;
}

std::unique_ptr<TH2D> Image::as_hist() const {
    if (header == nullptr) {throw except::invalid_operation("Error in Image::as_hist: Header must be initialized to use this method.");}
    std::unique_ptr<TH2D> hist = std::make_unique<TH2D>(utility::uid("hist").c_str(), "hist", header->nx, 0, header->cella_x, header->ny, 0, header->cella_y);
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
    for (unsigned int x = 0; x < N; x++) {
        bounds[x].min = 0;
        bounds[x].max = 0;
        bool min_set = false;
        for (unsigned int y = 0; y < M; y++) {
            if (index(x, y) < cutoff) {continue;}
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