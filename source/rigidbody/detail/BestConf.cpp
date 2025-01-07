/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/detail/BestConf.h>

using namespace ausaxs;

rigidbody::detail::BestConf::BestConf() = default;
rigidbody::detail::BestConf::BestConf(std::shared_ptr<grid::Grid> grid, std::vector<data::Water> waters, double chi2) noexcept : grid(std::move(grid)), waters(std::move(waters)), chi2(chi2) {}
rigidbody::detail::BestConf::~BestConf() = default;

