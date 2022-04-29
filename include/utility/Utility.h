#pragma once

#include <string>
#include <algorithm>

template<typename T>
/**
 * @brief Check if two numbers are approximately equal. 
 * 
 * @param v1 First value.
 * @param v2 Second value. 
 * @param abs Absolute tolerance. 
 * @param eps Relative tolerance. 
 */
bool approx(T v1, T v2, double abs = 1e-6, double eps = 0.01) {
    if (v1-abs > v2*(1+eps)) {return false;}
    if (v1+abs < v2*(1-eps)) {return false;}
    return true;
}

/**
 * @brief Remove spaces from both ends of a string. 
 *        Note that the input string is modified. 
 */
// std::string remove_spaces(std::string s) {
//     std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
//     s.erase(end_pos, s.end());
//     return s;
// }