#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Protein.h>

using namespace rigidbody;

OverlapConstraint::OverlapConstraint(Protein* protein) {
    this->protein = protein;
    initialize();
}

double OverlapConstraint::evaluate() const {
    auto distances = protein->get_histogram();
    double chi2 = 0;
    for (unsigned int i = 0; i < distances.size(); i++) {
        chi2 += std::pow(distances[i] - target[i], 2);
    }
    std::cout << "Overlap constraint: " << chi2 << std::endl;
    return chi2;
}

void OverlapConstraint::initialize() {
    auto initial_distances = protein->get_histogram();
    auto weight = [](double r) {
        return std::exp(-3*r);
    };

    // define the target distribution
    target = hist::Histogram(initial_distances.axis);
    auto axis = initial_distances.axis.as_vector();
    for (unsigned int i = 0; i < axis.size(); i++) {
        target[i] = weight(axis[i]) * initial_distances[i];
    }
}