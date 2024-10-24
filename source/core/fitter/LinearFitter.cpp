/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/LinearFitter.h>
#include <fitter/FitResult.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <math/CubicSpline.h>

using namespace ausaxs;
using namespace ausaxs::fitter;

LinearFitter::LinearFitter(const SimpleDataset& data) : data(data) {
    detail::LinearLeastSquares::data = data.y();
    detail::LinearLeastSquares::inv_sigma = data.yerr();
    for (unsigned i = 0; i < detail::LinearLeastSquares::inv_sigma.size(); ++i) {
        detail::LinearLeastSquares::inv_sigma[i] = 1./detail::LinearLeastSquares::inv_sigma[i];
    }
}

LinearFitter::LinearFitter(LinearFitter&&) noexcept = default;
LinearFitter& LinearFitter::operator=(LinearFitter&&) noexcept = default;

LinearFitter::LinearFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> model) : data(data), model(std::move(model)) {
    detail::LinearLeastSquares::data = data.y();
    detail::LinearLeastSquares::model = splice(this->model->debye_transform().get_counts());
    detail::LinearLeastSquares::inv_sigma = data.yerr();
    for (unsigned i = 0; i < detail::LinearLeastSquares::inv_sigma.size(); ++i) {
        detail::LinearLeastSquares::inv_sigma[i] = 1./detail::LinearLeastSquares::inv_sigma[i];
    }
}

void LinearFitter::refresh_model() {
    detail::LinearLeastSquares::model = splice(model->debye_transform().get_counts());
}

void LinearFitter::set_model(std::unique_ptr<hist::DistanceHistogram> model) {
    this->model = std::move(model);
    detail::LinearLeastSquares::model = splice(this->model->debye_transform().get_counts());
}

SimpleDataset LinearFitter::get_data() const {
    return data;
}

std::vector<double> LinearFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    math::CubicSpline s(model->get_q_axis(), ym);
    for (unsigned int i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x(i));
    }
    return Im;
}