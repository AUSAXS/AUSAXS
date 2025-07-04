// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Limit2D.h>

using namespace ausaxs;

Limit2D::Limit2D() noexcept = default;
Limit2D::Limit2D(const Limit& x, const Limit& y) noexcept : x(x), y(y) {}
Limit2D::Limit2D(double xmin, double xmax, double ymin, double ymax) noexcept : x(xmin, xmax), y(ymin, ymax) {}
bool Limit2D::empty() const noexcept {return x.empty() || y.empty();}