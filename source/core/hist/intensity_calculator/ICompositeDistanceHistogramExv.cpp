/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs;

Limit hist::ICompositeDistanceHistogramExv::get_excluded_volume_scaling_factor_limits() const {return {0.92, 1.08};}

Limit hist::ICompositeDistanceHistogramExv::get_solvent_density_scaling_factor_limits() const {return {0.5, 2};}

Limit hist::ICompositeDistanceHistogramExv::get_debye_waller_factor_limits() const {return {0.0, 5};}