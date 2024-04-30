/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Limit.h>

Limit::Limit() noexcept : min(0), max(0) {}

Limit::Limit(double min, double max) noexcept : min(min), max(max) {}

double Limit::span() const noexcept {return max-min;}

double Limit::center() const noexcept {return min + span()/2;}

void Limit::merge(const Limit& other) noexcept {
    min = std::min(min, other.min);
    max = std::max(max, other.max);
}

void Limit::expand(double percent) noexcept {
    double span = this->span();
    min -= span*percent;
    max += span*percent;
} 

Limit& Limit::operator-=(double c) noexcept {
    min-=c;
    max-=c;
    return *this;
}

Limit Limit::operator-(double c) const noexcept {return Limit(*this)-=c;}

Limit& Limit::operator+=(double c) noexcept {
    min+=c;
    max+=c;
    return *this;
}

Limit Limit::operator+(double c) const noexcept {return Limit(*this)+=c;}

bool Limit::operator==(const Limit& rhs) const noexcept {return min == rhs.min && max == rhs.max;}

bool Limit::operator!=(const Limit& rhs) const noexcept {return !operator==(rhs);}

std::string Limit::to_string(const std::string& prepend) const noexcept {return prepend + std::to_string(min) + " " + std::to_string(max);}

bool Limit::empty() const noexcept {return min == 0 && max == 0;}