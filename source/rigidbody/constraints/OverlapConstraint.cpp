/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <constants/ConstantsAxes.h>
#include <utility/Console.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

OverlapConstraint::OverlapConstraint(data::Molecule* protein) {
    this->protein = protein;
    initialize();
}

OverlapConstraint::~OverlapConstraint() = default;

void OverlapConstraint::set_overlap_function(std::function<double(double)> func) {
    overlap_function = std::move(func);
}

double OverlapConstraint::evaluate() const {
    if (target.empty()) [[unlikely]] {return 0;}
    auto current = protein->get_total_histogram()->get_total_counts();
    double chi2 = 0;
    for (unsigned int i = 1; i < target.size(); i++) { // skip the self-correlation bin
        chi2 += std::pow((current[i] - target[i])*weights[i], 2);
    }
    return chi2;
}

double OverlapConstraint::weight(double r) {
    return overlap_function(r);
}

void OverlapConstraint::initialize() {
    // define the target distribution
    auto hist = protein->get_histogram();
    target = hist->get_total_counts();
    axis = hist->get_d_axis();
    weights.resize(axis.size());

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