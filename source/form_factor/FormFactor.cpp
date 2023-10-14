#include <form_factor/FormFactor.h>
#include <data/record/Atom.h>

#include <cmath>

using namespace form_factor;

double FormFactor::evaluate(double q) const {
    double sum = 0;
    for (unsigned int i = 0; i < 5; ++i) {
        sum += a[i]*std::exp(-b[i]*q*q);
    }
    return (sum + c)/f0;
}