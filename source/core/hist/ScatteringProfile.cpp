// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/ScatteringProfile.h>

#include <dataset/Dataset.h>

#include <cassert>
#include <sstream>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringProfile::ScatteringProfile(std::vector<double> I, const Axis& q_axis) : I(std::move(I)), axis(q_axis) {
    assert(static_cast<int>(this->I.size()) == axis.bins && "ScatteringProfile: the intensity and q-axis must have the same size.");
}

const Axis& ScatteringProfile::get_axis() const {return axis;}

const std::vector<double>& ScatteringProfile::get_intensity() const {return I.data;}

int ScatteringProfile::size() const noexcept {return I.size();}

Dataset ScatteringProfile::as_dataset() const {
    return {axis.as_vector(), I.data};
}

std::string ScatteringProfile::to_string() const {
    std::stringstream ss;
    auto q = axis.as_vector();
    for (int i = 0; i < size(); ++i) {
        ss << q[i] << " " << I[i] << std::endl;
    }
    return ss.str();
}

ScatteringProfile& ScatteringProfile::operator+=(const ScatteringProfile& rhs) {I += rhs.I; return *this;}
ScatteringProfile& ScatteringProfile::operator-=(const ScatteringProfile& rhs) {I -= rhs.I; return *this;}
ScatteringProfile& ScatteringProfile::operator*=(double rhs) {I *= rhs; return *this;}
double& ScatteringProfile::operator[](int i) {return I[i];}
double ScatteringProfile::operator[](int i) const {return I[i];}
bool ScatteringProfile::operator==(const ScatteringProfile& rhs) const {return I == rhs.I && axis == rhs.axis;}

ScatteringProfile hist::operator+(const ScatteringProfile& lhs, const ScatteringProfile& rhs) {
    ScatteringProfile result(lhs);
    result += rhs;
    return result;
}

ScatteringProfile hist::operator-(const ScatteringProfile& lhs, const ScatteringProfile& rhs) {
    ScatteringProfile result(lhs);
    result -= rhs;
    return result;
}

ScatteringProfile hist::operator*(const ScatteringProfile& lhs, double rhs) {
    ScatteringProfile result(lhs);
    result *= rhs;
    return result;
}
