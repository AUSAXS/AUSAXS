#include <hist/detail/MasterHistogram.h>

using namespace hist::detail;

MasterHistogram::MasterHistogram(const std::vector<double>& p_base, const Axis& axis) : Histogram(p_base, axis), base(std::move(p_base)) {}

MasterHistogram& MasterHistogram::operator+=(const PartialHistogram& rhs) {
    p += Vector(rhs.get_counts());
    return *this;
}

MasterHistogram& MasterHistogram::operator-=(const PartialHistogram& rhs) {
    p -= Vector(rhs.get_counts());
    return *this;
}