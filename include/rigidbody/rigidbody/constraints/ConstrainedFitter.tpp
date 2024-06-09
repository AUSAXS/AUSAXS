#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/Exceptions.h>

using namespace fitter;

template<fitter::fitter_t T>
ConstrainedFitter<T>::ConstrainedFitter(ConstrainedFitter<T>&& other) : T(std::move(other)), constraints(std::move(other.constraints)) {}

template<fitter::fitter_t T>
double ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    return T::chi2(params) + constraints->evaluate();
}

template<fitter::fitter_t T>
void ConstrainedFitter<T>::set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints) {
    this->constraints = constraints;
}

template<fitter::fitter_t T>
observer_ptr<rigidbody::constraints::ConstraintManager> ConstrainedFitter<T>::get_constraint_manager() {
    return this->constraints.get();
}