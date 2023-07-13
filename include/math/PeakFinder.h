/**
MIT License

Copyright (c) 2019 Clayder Gonzalez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
**/

#pragma once

#include <vector>

namespace peak_finder {
    constexpr double eps = 1e-16;

   /**
    * @brief Find peaks in a signal.
    * 
    * @param x0 Input data to find peaks in.
    * @param includeEndpoints If true the endpoints will be included as possible extrema. Otherwise they will not be included.
    * @param extrema 1 if maxima are desired, -1 if minima are desired.
    * 
    * @return Indices of peaks in the data.
    */
    std::vector<unsigned int> find_peaks(std::vector<double> data, bool includeEndpoints = true, double extrema = 1);
}