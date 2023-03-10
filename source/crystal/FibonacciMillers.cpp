#include <crystal/miller/FibonacciMillers.h>

#include <math.h>

using namespace crystal;

FibonacciMillers::FibonacciMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

std::vector<Miller> FibonacciMillers::generate() const {
    std::vector<Miller> bases = pick_directions(generate_independent_bases());
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

#include <iostream>
#include <fstream>
std::vector<Miller> FibonacciMillers::pick_directions(const std::vector<Miller>& basis) const {
    std::vector<Vector3<double>> normed(basis.size());
    std::transform(basis.begin(), basis.end(), normed.begin(), [](const auto& miller) {return miller.normalize();});

    // calculate smallest distance between miller indices
    double res = 1e10;
    for (const auto& miller : normed) {
        for (const auto& other : normed) {
            if (miller == other) {continue;}
            double distance = miller.distance(other);
            if (distance < res) {
                res = distance;
            }
        }
    }

    // sqrt(n)*res = 1
    // n = 1/res^2

    // generate a fibonacci sphere with roughly the same resolution as the smallest distance
    auto fib = generate_fibonacci_sphere(1./(res*res));

    // find the closest miller index for each point on the sphere
    std::vector<Miller> directions;
    std::vector<Vector3<double>> v;
    for (const auto& point : fib) {
        double min_distance = 1e10;
        unsigned int index = 0;
        for (unsigned int i = 0; i < basis.size(); i++) {
            auto norm = normed[i];
            double distance = std::sqrt(std::pow(point.x() - norm.x(), 2) + std::pow(point.y() - norm.y(), 2) + std::pow(point.z() - norm.z(), 2));
            if (distance < min_distance) {
                min_distance = distance;
                index = i;
            }
        }
        if (min_distance < res) {
            directions.push_back(basis[index]);
            v.push_back(normed[index]);
        }
    }

    {
        std::ofstream file("fibonacci_sphere.py");
        file << "import numpy as np" << std::endl;
        file << "import matplotlib.pyplot as plt" << std::endl;
        file << "fig = plt.figure()" << std::endl;
        file << "x = np.array([";
        for (const auto& point : v) {
            file << point.x() << ", ";
        }
        file << "])" << std::endl;
        file << "y = np.array([";
        for (const auto& point : v) {
            file << point.y() << ", ";
        }
        file << "])" << std::endl;
        file << "z = np.array([";
        for (const auto& point : v) {
            file << point.z() << ", ";
        }
        file << "])" << std::endl;
        file << "ax = fig.add_subplot(111, projection='3d')" << std::endl;
        file << "ax.scatter(x, y, z)" << std::endl;
        file << "plt.show()" << std::endl;
        file.close();
    }

    return directions;
}

int FibonacciMillers::estimate_n(double resolution) const {
    // return resolution/(2*M_PI)*phi;
    return 1000;
}

std::vector<Vector3<double>> FibonacciMillers::generate_fibonacci_sphere(int n) const {
    std::cout << "Generating fibonacci sphere with " << n << " points" << std::endl;
    std::vector<Vector3<double>> points(n);

    for (int i = 0; i < n; i++) {
        double theta = 2*M_PI*i/phi;
        double z = 1 - 2*(i + 0.5)/n;
        double r = sqrt(1 - z*z);
        points[i] = Vector3<double>(r*cos(theta), r*sin(theta), z);
    }

    return points;
}
