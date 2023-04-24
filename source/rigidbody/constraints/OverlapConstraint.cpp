#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Protein.h>
#include <settings/GeneralSettings.h>

using namespace rigidbody;

OverlapConstraint::OverlapConstraint(Protein* protein) {
    this->protein = protein;
    initialize();
}

double OverlapConstraint::evaluate() const {
    auto distances = protein->get_histogram();
    double chi2 = 0;
    for (unsigned int i = 0; i < target.size(); i++) {
        chi2 += std::pow((distances[i] - target[i])*weights[i], 2);
    }
    // std::cout << "Overlap constraint: " << chi2 << std::endl;
    return chi2;
}

void OverlapConstraint::initialize() {
    auto weight = [](double r) {
        return std::exp(-5*r);
    };

    // define the target distribution
    target = protein->get_total_histogram();
    auto axis = target.axis.as_vector();
    weights = hist::Histogram(axis);

    // calculate the weights & determine how much of the histogram we have to use
    for (unsigned int i = 0; i < axis.size(); ++i) {
        weights[i] = weight(axis[i]);
        if (weights[i] < 1e-3) {weights[i] = 0;}
    }

    unsigned int i = axis.size()-1;
    for (; i > 0; i--) {
        if (weights[i] != 0) {break;}
    }

    target.resize(i);
    weights.resize(i);

    if (settings::general::verbose) {
        std::cout << "\tOverlap constraint initialized. The distance range [0, " << axis[i] << "]Å will be used for calculating the overlap penalty." << std::endl;
    }
}