#include <rigidbody/ConstrainedFitter.h>

using namespace fitter;

double ConstrainedFitter::chi2(const std::vector<double>& params) {
    double chi2 = HydrationFitter::chi2(params);
    for (const rigidbody::Constraint& constraint : constraints) {
        chi2 += constraint.evaluate();
    }
    return chi2;
}

void ConstrainedFitter::add_constraint(const rigidbody::Constraint& constraint) {
    constraints.push_back(constraint);
}

void ConstrainedFitter::add_constraint(rigidbody::Constraint&& constraint) {
    constraints.push_back(std::move(constraint));
}

void ConstrainedFitter::set_constraints(const std::vector<rigidbody::Constraint>& constraints) {
    this->constraints = constraints;
}

void ConstrainedFitter::set_constraints(std::vector<rigidbody::Constraint>&& constraints) {
    this->constraints = std::move(constraints);
}

const std::vector<rigidbody::Constraint>& ConstrainedFitter::get_constraints() const {
    return constraints;
}