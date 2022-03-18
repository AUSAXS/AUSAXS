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