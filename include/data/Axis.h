#pragma once

#include <vector>
#include <array>
#include <initializer_list>
#include <math.h>

using std::vector;

struct Limit {
    Limit(double min, double max) : min(min), max(max) {}

    double span() const {return max-min;}

    double min, max;
};

class Limit3D {
  public:
    Limit3D(const Limit& x, const Limit& y, const Limit& z) : x(x), y(y), z(z) {}
    Limit3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) : x(xmin, xmax), y(ymin, ymax), z(zmin, zmax) {}

    Limit x, y, z;  
};

class Axis {
  public:
    Axis() : bins(0), min(0), max(0) {}
    Axis(const Limit& limits, int bins) : bins(bins), min(limits.min), max(limits.max) {}
    Axis(const int bins, const double xmin, const double xmax) : bins(bins), min(xmin), max(xmax)  {}

    Axis& operator=(std::initializer_list<double> list) {
        std::vector<double> d = list;
        bins = std::round(d[0]); 
        min = d[1];
        max = d[2];
        return *this;
    }

    std::string to_string() const {
        return "Axis: (" + std::to_string(min) + ", " + std::to_string(max) + ") with " + std::to_string(bins) + " bins";
    }

    bool operator==(const Axis& rhs) const {
        if (bins != rhs.bins) {return false;}
        if (min != rhs.min) {return false;}
        if (max != rhs.max) {return false;}
        return true;
    }

    bool operator!=(const Axis& rhs) const {return !operator==(rhs);}

    friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {os << axis.to_string(); return os;}

    double width() const {
        return (max-min)/bins;
    }

    int bins;
    double min, max;
};

class Axis3D {
  public:
    Axis3D() {}
    Axis3D(const Axis3D& axis) : x(axis.x), y(axis.y), z(axis.z) {}
    Axis3D(const Axis& x, const Axis& y, const Axis& z) : x(x), y(y), z(z) {}
    Axis3D(const Limit3D& limits, double width) : x(limits.x, limits.x.span()/width), y(limits.y, limits.y.span()/width), z(limits.z, limits.z.span()/width) {}
    Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int bins) : x(bins, xmin, xmax), y(bins, ymin, ymax), z(bins, zmin, zmax) {}
    Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double width) : x((xmax-xmin)/width, xmin, xmax), y((ymax-ymin)/width, ymin, ymax), z((zmax-zmin)/width, zmin, zmax) {}
    Axis3D(const vector<int>& min, const vector<int>& max, double width) : x((max[0]-min[0])/width, min[0], max[0]), y((max[1]-min[1])/width, min[1], max[1]), z((max[2]-min[2])/width, min[2], max[2]) {}

    Axis3D& operator=(const Axis3D& rhs) {
      x = rhs.x;
      y = rhs.y;
      z = rhs.z;
      return *this;
    }

    bool operator==(const Axis3D& rhs) const {
      if (x != rhs.x) {return false;}
      if (y != rhs.y) {return false;}
      if (z != rhs.z) {return false;}
      return true;
    }

    bool operator!=(const Axis3D& rhs) const {return !operator==(rhs);}

    Axis x, y, z;
};