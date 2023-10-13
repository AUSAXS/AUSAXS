#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Molecule.h>
#include <hist/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>

using namespace rigidbody;

OverlapConstraint::OverlapConstraint() = default;

OverlapConstraint::OverlapConstraint(data::Molecule* protein) {
    this->protein = protein;
    initialize();
}

OverlapConstraint::~OverlapConstraint() = default;

double OverlapConstraint::evaluate() const {
    if (target.size() == 0) [[unlikely]] {return 0;}
    auto distances = protein->get_histogram()->get_d_axis();
    double chi2 = 0;
    for (unsigned int i = 0; i < target.size(); i++) {
        chi2 += std::pow((distances[i] - target[i])*weights[i], 2);
    }
    // std::cout << "Overlap constraint: " << chi2 << std::endl;
    return chi2;
}

double OverlapConstraint::weight(double r) {
    return std::exp(-5*r);
}

void OverlapConstraint::initialize() {
    // define the target distribution
    target = *protein->get_total_histogram().get();
    auto axis = target.get_axis().as_vector();
    weights = hist::Histogram(axis);

    // calculate the weights and reduce their precision
    for (unsigned int i = 0; i < axis.size(); ++i) {
        weights[i] = weight(axis[i]);
        if (weights[i] < 1e-3) {weights[i] = 0;}
    }

    // find the last non-zero weight
    unsigned int i = axis.size()-1;
    for (; i > 0; i--) {
        if (weights[i] != 0) {break;}
    }

    // resize the histograms to the last non-zero weight
    target.resize(i);
    weights.resize(i);

    if (settings::general::verbose) {
        std::cout << "\tOverlap constraint initialized. The distance range [0, " << axis[i] << "]Å will be used for calculating the overlap penalty." << std::endl;
        if (target.size() < 5) {console::print_warning("\tWarning: Only " + std::to_string(target.size()) + " bins will be used for calculating the overlap penalty. Consider decreasing the bin size in the histogram settings.");}
    }
}

bool OverlapConstraint::operator==(const OverlapConstraint& other) const = default;