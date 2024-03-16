/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

// #include <hist/ScatteringHistogramOld.h>
// #include <settings/HistogramSettings.h>
// #include <utility/Constants.h>
// #include <utility/Axis.h>

// #include <utility>
// #include <math.h>
// #include <algorithm>
// #include <iostream>

// using namespace hist;

// ScatteringHistogram::ScatteringHistogram() = default;

// ScatteringHistogram::ScatteringHistogram(const ScatteringHistogram&& sh) noexcept : Histogram(sh.p, sh.axis), p_pp(sh.p_pp), p_hh(sh.p_hh), p_hp(sh.p_hp) {
//     setup();
// }

// ScatteringHistogram::ScatteringHistogram(const ScatteringHistogram& sh) : Histogram(sh.p, sh.axis), p_pp(sh.p_pp), p_hh(sh.p_hh), p_hp(sh.p_hp) {
//     setup();
// }

// ScatteringHistogram::ScatteringHistogram(const std::vector<double>& p_pp, const std::vector<double>& p_hh, const std::vector<double>& p_hp, const std::vector<double>& p_tot, const Axis& axis)
//     : Histogram(p_tot, axis), p_pp(p_pp, axis), p_hh(p_hh, axis), p_hp(p_hp, axis) {setup();}

// ScatteringHistogram::~ScatteringHistogram() = default;

// void ScatteringHistogram::setup() {
//     Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
//     q = debye_axis.as_vector();
//     d = axis.as_vector(true);
//     d[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms

//     sinqd_table = std::make_unique<table::DebyeLookupTable>(q, d);
// }

// void ScatteringHistogram::apply_water_scaling_factor(const double& k) {
//     double k2 = std::pow(k, 2);
//     for (unsigned int i = 0; i < axis.bins; ++i) {p[i] = p_pp[i] + k*p_hp[i] + k2*p_hh[i];} // p = p_tot, inherited from Histogram
// }

// Histogram ScatteringHistogram::plot_debye_scattering() const {
//     SimpleDataset data = calc_debye_scattering_intensity();
//     return Histogram(data.y(), Axis(q[0], q.back(), q.size()));
// }

// SimpleDataset ScatteringHistogram::calc_debye_scattering_intensity(const std::vector<double>& q) const {
//     // calculate the scattering intensity based on the Debye equation
//     std::vector<double> Iq(q.size(), 0);
//     for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
//         for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
//             double qd = q[i]*d[j];
//             if (qd < 1e-6) {Iq[i] += 1;}
//             else {Iq[i] += p[j]*std::sin(qd)/qd;}
//         }
//         Iq[i] *= std::exp(-q[i]*q[i]); // form factor
//     }
//     return SimpleDataset(q, Iq, "q", "I");
// }

// SimpleDataset ScatteringHistogram::calc_debye_scattering_intensity() const {
//     // calculate the Debye scattering intensity
//     Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

//     // calculate the scattering intensity based on the Debye equation
//     std::vector<double> Iq(debye_axis.bins, 0);
//     for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
//         for (unsigned int j = 0; j < p.size(); j++) { // iterate through the distance histogram
//             Iq[i] += p[j]*sinqd_table->lookup(i, j);
//         }
//         Iq[i] *= std::exp(-q[i]*q[i]); // form factor
//     }
//     return SimpleDataset(q, Iq, "q", "I");
// }

// Histogram ScatteringHistogram::plot_guinier_approx() const {
//     std::vector<double> Iq = calc_guinier_approx().col("logI");
//     return Histogram(Iq, Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins));
// }

// double ScatteringHistogram::calc_guinier_gyration_ratio_squared() const {
//     double num = 0, denom = 0;
//     for (unsigned int i = 0; i < p.size(); ++i) {
//         num += p[i]*std::pow(d[i], 2);
//         denom += 2*p[i];
//     }
//     return num/denom;
// }

// SimpleDataset ScatteringHistogram::calc_guinier_approx() const {
//     double Rg2 = calc_guinier_gyration_ratio_squared();

//     Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
//     std::vector<double> Iq(debye_axis.bins, 0);
//     for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
//         Iq[i] = std::exp(-std::pow(q[i], 2)*Rg2/3);
//     }

//     return SimpleDataset(q, Iq, "q", "I");
// }

// std::string ScatteringHistogram::to_string() const noexcept {
//     return calc_debye_scattering_intensity().to_string();
// }

// ScatteringHistogram& ScatteringHistogram::operator=(const ScatteringHistogram& h) {
//     *this = std::move(ScatteringHistogram(h));
//     return *this;
// }

// ScatteringHistogram& ScatteringHistogram::operator=(ScatteringHistogram&& h) = default;

// bool ScatteringHistogram::operator==(const ScatteringHistogram& h) const {
//     return Histogram::operator==(h) && p_pp == h.p_pp && p_hh == h.p_hh && p_hp == h.p_hp;
// }

// ScatteringHistogram& ScatteringHistogram::operator+=(const ScatteringHistogram& rhs) {
//     Histogram::operator+=(rhs);
//     p_pp += rhs.p_pp;
//     p_hh += rhs.p_hh;
//     p_hp += rhs.p_hp;
//     return *this;
// }

// ScatteringHistogram& ScatteringHistogram::operator-=(const ScatteringHistogram& rhs) {
//     Histogram::operator-=(rhs);
//     p_pp -= rhs.p_pp;
//     p_hh -= rhs.p_hh;
//     p_hp -= rhs.p_hp;
//     return *this;
// }

// ScatteringHistogram& ScatteringHistogram::operator*=(double rhs) {
//     Histogram::operator*=(rhs);
//     p_pp *= rhs;
//     p_hh *= rhs;
//     p_hp *= rhs;
//     return *this;
// }

// ScatteringHistogram hist::operator+(const ScatteringHistogram& lhs, const ScatteringHistogram& rhs) {
//     ScatteringHistogram result(lhs);
//     result += rhs;
//     return result;
// }

// ScatteringHistogram hist::operator-(const ScatteringHistogram& lhs, const ScatteringHistogram& rhs) {
//     ScatteringHistogram result(lhs);
//     result -= rhs;
//     return result;
// }

// ScatteringHistogram hist::operator*(const ScatteringHistogram& lhs, double rhs) {
//     ScatteringHistogram result(lhs);
//     result *= rhs;
//     return result;
// }

// // void ScatteringHistogram::extend_axis(double qmax) {
// //     Histogram::extend_axis(qmax);
// //     p_pp.extend_axis(qmax);
// //     p_hh.extend_axis(qmax);
// //     p_hp.extend_axis(qmax);
// // }