#include <rigidbody/ConstrainedFitter.h>

using namespace fitter;

template<fitter_t T>
double ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    double chi2 = T::chi2(params);
    for (const auto& constraint : constraints) {
        chi2 += constraint->evaluate();
    }
    return chi2;
}

template<fitter_t T>
void ConstrainedFitter<T>::add_constraint(std::shared_ptr<rigidbody::Constraint> constraint) {
    constraints.push_back(std::move(constraint));
}

template<fitter_t T>
void ConstrainedFitter<T>::add_constraint(rigidbody::Constraint&& constraint) {
    constraints.push_back(std::make_shared<rigidbody::Constraint>(std::move(constraint)));
}

template<fitter_t T>
void ConstrainedFitter<T>::set_constraints(std::vector<std::shared_ptr<rigidbody::Constraint>> constraints) {
    this->constraints = constraints;
}

template<fitter_t T>
const std::vector<std::shared_ptr<rigidbody::Constraint>>& ConstrainedFitter<T>::get_constraints() const {
    return constraints;
}