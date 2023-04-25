#pragma once

#include <math/Vector3.h>
#include <utility/Limit.h>

#include <vector>
#include <string>

struct CoordinateSystem {
    CoordinateSystem() = default;
    CoordinateSystem(const Vector3<float>& origin, const Vector3<float>& x, const Vector3<float>& y, const Vector3<float>& z) : origin(origin), x(x), y(y), z(z) {}

    std::string to_string() const {
        std::string result = "origin: " + origin.to_string() + "\n";
        result += "x: " + x.to_string() + "\n";
        result += "y: " + y.to_string() + "\n";
        result += "z: " + z.to_string() + "\n";
        return result;
    }

    Vector3<float> origin, x, y, z;    
};

namespace em {
    class Box {
        public: 
            Box() = default;
            Box(Limit x, Limit y, Limit z) : x(x), y(y), z(z) {}

            void save(std::string filename) const;

            void distances() const;

            /**
             * @brief Generate a coordinate system spanning the @a contents of this box.
             * 
             * The axes are defined with respect to origo. 
             */
            CoordinateSystem spanning_coordinate_system() const;

            Limit x, y, z;
            std::vector<Vector3<float>> atoms;
    };
}