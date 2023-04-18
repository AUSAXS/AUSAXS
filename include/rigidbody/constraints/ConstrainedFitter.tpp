#include <rigidbody/constraints/ConstrainedFitter.h>

using namespace fitter;

template<fitter_t T>
double ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    double chi2 = T::chi2(params);
    for (const auto& constraint : distance_constraints) {
        chi2 += constraint->evaluate();
    }
    chi2 += overlap_constraint->evaluate();
    return chi2;
}

template<fitter_t T>
void ConstrainedFitter<T>::add_constraint(std::shared_ptr<rigidbody::DistanceConstraint> constraints) {
    distance_constraints.push_back(std::move(constraints));
}

template<fitter_t T>
void ConstrainedFitter<T>::add_constraint(rigidbody::DistanceConstraint&& constraint) {
    constraints.push_back(std::make_shared<rigidbody::DistanceConstraint>(std::move(constraint)));
}

template<fitter_t T>
void ConstrainedFitter<T>::set_constraints(std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> constraints) {
    this->constraints = constraints;
}

template<fitter_t T>
void ConstrainedFitter<T>::set_constraints(std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> constraints) {
    this->constraints = constraints;
}

template<fitter_t T>
const std::vector<std::shared_ptr<rigidbody::DistanceConstraint>>& ConstrainedFitter<T>::get_constraints() const {
    return distance_constraints;
}