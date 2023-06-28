#include <utility/Limit3D.h>

Limit3D::Limit3D() noexcept = default;

Limit3D::Limit3D(const Limit& x, const Limit& y, const Limit& z) noexcept : x(x), y(y), z(z) {}

Limit3D::Limit3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) noexcept : x(xmin, xmax), y(ymin, ymax), z(zmin, zmax) {}

bool Limit3D::empty() const noexcept {return x.empty() || y.empty() || z.empty();}