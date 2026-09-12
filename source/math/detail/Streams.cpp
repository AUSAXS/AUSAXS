// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/Vector.h>
#include <math/Vector3.h>
#include <math/slices/Slice.h>

#include <ostream>

namespace ausaxs {
    template<numeric T>
    std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {os << v.to_string(); return os;}

    template<numeric T>
    std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {os << v.to_string(); return os;}

    template<numeric T, container_type Container>
    std::ostream& operator<<(std::ostream& os, const Slice<T, Container>& v) {os << v.to_string(); return os;}

    template std::ostream& operator<< <double>(std::ostream&, const Vector<double>&);
    template std::ostream& operator<< <int>   (std::ostream&, const Vector<int>&);

    template std::ostream& operator<< <double>(std::ostream&, const Vector3<double>&);
    template std::ostream& operator<< <int>   (std::ostream&, const Vector3<int>&);
    template std::ostream& operator<< <float> (std::ostream&, const Vector3<float>&);

    template std::ostream& operator<< <double, std::vector<double>&>      (std::ostream&, const Slice<double, std::vector<double>&>&);
    template std::ostream& operator<< <double, const std::vector<double>&>(std::ostream&, const Slice<double, const std::vector<double>&>&);
    template std::ostream& operator<< <float,  std::vector<float>&>       (std::ostream&, const Slice<float, std::vector<float>&>&);
    template std::ostream& operator<< <float,  const std::vector<float>&> (std::ostream&, const Slice<float, const std::vector<float>&>&);
}
