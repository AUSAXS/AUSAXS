/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Limit2D.h>

Limit2D::Limit2D() noexcept = default;
Limit2D::Limit2D(const Limit& x, const Limit& y) noexcept : x(x), y(y) {}
Limit2D::Limit2D(double xmin, double xmax, double ymin, double ymax) noexcept : x(xmin, xmax), y(ymin, ymax) {}
bool Limit2D::empty() const noexcept {return x.empty() || y.empty();}