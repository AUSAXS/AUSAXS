// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/detail/EMFitResult.h>
#include <utility/Utility.h>

using namespace ausaxs::fitter;

std::string EMFitResult::to_string() const noexcept {
    std::stringstream ss;
    ss << FitResult::to_string();

                    ss << "\n| Cutoff corresponds to PyMOL level " << utility::print_element(level, 12) << "           |";
    if (mass != 0) {ss << "\n|                  and to a mass of " << utility::print_element(mass, 12)  << " kDa       |";}
                    ss << "\n+-----------------------------------"                  "------------"                  "-----------+";

    return ss.str();
}