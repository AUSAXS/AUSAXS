#include <rigidbody/constraints/ConstrainedFitter.h>

template<fitter::fitter_t T>
double fitter::ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    std::cout << "ConstrainedFitter::chi2() called." << std::endl;
    return T::chi2(params) + constraints->evaluate();
}

template<fitter::fitter_t T>
void fitter::ConstrainedFitter<T>::set_constraint_manager(std::shared_ptr<rigidbody::ConstraintManager> constraints) {
    this->constraints = constraints;
}