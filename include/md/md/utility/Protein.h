#pragma once

#include <string>
#include <vector>

namespace ausaxs::md {
    struct Vector3 {
        double x, y, z;
    };

    class Protein {
        public:
            Protein(const std::string& filename) {
                atoms = parse(filename);
            }

            double maximum_distance() const;

        private:
            std::vector<Vector3> atoms;
            std::vector<Vector3> parse(const std::string& filename);
    };
}