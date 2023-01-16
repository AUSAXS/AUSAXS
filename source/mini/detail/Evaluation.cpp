#include <mini/detail/Evaluation.h>

#include <string>

std::string mini::Evaluation::to_string() const {
    std::string s;
    for (auto val : vals) {
        s += std::to_string(val) + " ";
    }
    s += std::to_string(fval);
    return s;
}