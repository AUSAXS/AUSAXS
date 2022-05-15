#include <minimizer/Minimizer.h>

#include <functional>

std::vector<double> Minimizer::Result::parameters() const {}

std::vector<double> Minimizer::Result::errors() const {}

void Minimizer::set_function(double(&f)(double*)) {
    function = std::bind(f, std::placeholders::_1);
}