/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/EMGrid.h>
#include <settings/GridSettings.h>

using namespace ausaxs::em::grid;

double EMGrid::get_atomic_radius(form_factor::form_factor_t) const {
    return settings::grid::min_exv_radius;
}