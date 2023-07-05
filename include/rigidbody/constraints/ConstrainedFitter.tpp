#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/Exceptions.h>
#include <Symbols.h>

template<fitter::fitter_t T>
double fitter::ConstrainedFitter<T>::chi2(const std::vector<double>& params) {
    #ifdef DEBUG
        if (constraints == nullptr) {throw except::nullptr_error("ConstrainedFitter::chi2: Constraint manager is not set.");}
    #endif

    return T::chi2(params) + constraints->evaluate();
}

template<fitter::fitter_t T>
void fitter::ConstrainedFitter<T>::set_constraint_manager(std::shared_ptr<rigidbody::ConstraintManager> constraints) {
    this->constraints = constraints;
}

template<fitter::fitter_t T>
rigidbody::ConstraintManager* fitter::ConstrainedFitter<T>::get_constraint_manager() {
    #ifdef DEBUG
        if (constraints == nullptr) {throw except::nullptr_error("ConstrainedFitter::get_constraint_manager: Constraint manager is not set.");}
    #endif

    return this->constraints.get();
}