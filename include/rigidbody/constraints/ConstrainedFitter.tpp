#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>

template<fitter::fitter_t T>
double fitter::ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    return T::chi2(params) + constraints->evaluate();
}

template<fitter::fitter_t T>
void fitter::ConstrainedFitter<T>::set_constraint_manager(std::shared_ptr<rigidbody::ConstraintManager> constraints) {
    this->constraints = constraints;
}

template<fitter::fitter_t T>
std::shared_ptr<rigidbody::ConstraintManager> fitter::ConstrainedFitter<T>::get_constraint_manager() {
    return constraints;
} 