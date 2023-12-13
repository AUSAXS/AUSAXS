#include <crystal/miller/FibonacciMillers.h>
#include <crystal/miller/Miller.h>
#include <crystal/Fval.h>
#include <settings/CrystalSettings.h>
#include <constants/ConstantsMath.h>

#include <math.h>
#include <iostream>
#include <fstream>

using namespace crystal;

FibonacciMillers::FibonacciMillers(unsigned int h, unsigned int k, unsigned int l) : ReducedMillers(h, k, l) {}

std::vector<Miller> FibonacciMillers::generate() const {
    std::vector<Miller> bases = pick_directions(generate_independent_bases(settings::crystal::reduced::basis_q));
    std::vector<Miller> millers;

    // now generate all millers indices
    // we can do this by multiplying the base pairs with integers
    for (const auto& base : bases) {
        int multiplier = 1;
        while (multiplier++ < 10000) { // hard limit to prevent infinite loop
            Miller miller(base.h*multiplier, base.k*multiplier, base.l*multiplier);
            double q = crystal::Fval::Q(miller).norm();
            if (q > settings::crystal::max_q) {break;}
            millers.emplace_back(miller);
        }
    }

    return millers;
}

std::vector<Miller> FibonacciMillers::pick_directions(const std::vector<Miller>& basis) const {
    auto copy = basis;
    return pick_directions(std::move(copy));
}

std::vector<Miller> FibonacciMillers::pick_directions(std::vector<Miller>&& bases) const {
    std::vector<Vector3<double>> v = generate_fibonacci_sphere(200);
    std::vector<Vector3<double>> normed(bases.size());
    std::transform(bases.begin(), bases.end(), normed.begin(), [](const Miller& m) {return m.normalize();});

    std::vector<Miller> directions;
    std::vector<Vector3<double>> directions2;
    for (const auto& point : v) {
        double min_distance = 1000000;
        unsigned int closest = 0;
        for (unsigned int i = 0; i < normed.size(); i++) {
            double distance = point.distance2(normed[i]);
            if (distance < min_distance) {
                min_distance = distance;
                closest = i;
            }
        }
        if (min_distance < 0.1) {
            directions.push_back(bases[closest]);
            directions2.push_back(normed[closest]);
        }
    }

    std::cout << "Found " << directions.size() << " directions. Expected " << v.size() << "." << std::endl;

    // {
    //     std::ofstream file1("fs1.py");
    //     file1 << "import numpy as np" << std::endl;
    //     file1 << "import matplotlib.pyplot as plt" << std::endl;
    //     file1 << "fig = plt.figure()" << std::endl;
    //     file1 << "x = np.array([";
    //     for (const auto& point : directions) {
    //         file1 << point.h << ", ";
    //     }
    //     file1 << "])" << std::endl;
    //     file1 << "y = np.array([";
    //     for (const auto& point : directions) {
    //         file1 << point.k << ", ";
    //     }
    //     file1 << "])" << std::endl;
    //     file1 << "z = np.array([";
    //     for (const auto& point : directions) {
    //         file1 << point.l << ", ";
    //     }
    //     file1 << "])" << std::endl;
    //     file1 << "ax = fig.add_subplot(111, projection='3d')" << std::endl;
    //     file1 << "ax.scatter(x, y, z)" << std::endl;
    //     file1 << "plt.show()" << std::endl;
    //     file1.close();

    //     std::ofstream file2("fs2.py");
    //     file2 << "import numpy as np" << std::endl;
    //     file2 << "import matplotlib.pyplot as plt" << std::endl;
    //     file2 << "fig = plt.figure()" << std::endl;
    //     file2 << "x = np.array([";
    //     for (const auto& point : directions2) {
    //         file2 << point.x() << ", ";
    //     }
    //     file2 << "])" << std::endl;
    //     file2 << "y = np.array([";
    //     for (const auto& point : directions2) {
    //         file2 << point.y() << ", ";
    //     }
    //     file2 << "])" << std::endl;
    //     file2 << "z = np.array([";
    //     for (const auto& point : directions2) {
    //         file2 << point.z() << ", ";
    //     }
    //     file2 << "])" << std::endl;
    //     file2 << "ax = fig.add_subplot(111, projection='3d')" << std::endl;
    //     file2 << "ax.scatter(x, y, z)" << std::endl;
    //     file2 << "plt.show()" << std::endl;
    //     file2.close();
    // }

    return directions;
}

int FibonacciMillers::estimate_n(double) const {
    throw std::runtime_error("Not implemented");
    // return resolution/(2*constants::pi)*phi;
    // return 1000;
}

std::vector<Vector3<double>> FibonacciMillers::generate_fibonacci_sphere(int n) const {
    std::cout << "Generating fibonacci sphere with " << n << " points" << std::endl;
    std::vector<Vector3<double>> points(n);

    for (int i = 0; i < n; i++) {
        double theta = 2*constants::pi*i/phi;
        double z = 1 - 2*(i + 0.5)/n;
        double r = sqrt(1 - z*z);
        points[i] = Vector3<double>(r*cos(theta), r*sin(theta), z);
    }

    return points;
}
