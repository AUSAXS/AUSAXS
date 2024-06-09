/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <math/MovingAverager.h>

#include <stdexcept>

void MovingAverage::validate_input(unsigned int N, unsigned int window_size) {
    if (N < window_size) {
        throw std::invalid_argument("MovingAverager::average: Window size is larger than data size.");
    }

    if (window_size % 2 == 0) {
        throw std::invalid_argument("Error in MovingAverager::validate_input: Window_size must be odd");
    }
}
