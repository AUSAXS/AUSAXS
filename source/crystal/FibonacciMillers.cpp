#include <crystal/miller/FibonacciMillers.h>

#include <math.h>
#include <iostream>

using namespace crystal;

FibonacciMillers::FibonacciMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

std::vector<Miller> FibonacciMillers::generate() const {
    std::vector<Miller> bases = generate_independent_bases();
    std::vector<Miller> millers;

    // now generate all millers indices
    // we can do this by multiplying the base pairs with integers
    for (const auto& base : bases) {
        int multiplier = 0;
        while (multiplier++ < 100000) { // hard limit to prevent infinite loop
            if (base.h*multiplier > h || base.k*multiplier > k || base.l*multiplier > l) {break;}
            millers.emplace_back(base.h*multiplier, base.k*multiplier, base.l*multiplier);
        }
    }

    return millers;
}

std::vector<Vector3<double>> FibonacciMillers::generate_fibonacci_sphere(int n) const {
    std::vector<Vector3<double>> points(n);

    double phi = (1 + sqrt(5))/2;
    for (int i = 0; i < n; i++) {
        double theta = 2*M_PI*i/phi;
        double z = 1 - 2*(i + 0.5)/n;
        double r = sqrt(1 - z*z);
        points[i] = Vector3<double>(r*cos(theta), r*sin(theta), z);
    }

    return points;
}
