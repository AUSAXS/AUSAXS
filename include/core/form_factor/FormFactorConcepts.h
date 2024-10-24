#include <concepts>

namespace ausaxs {
    template<typename T>
        concept FormFactorType = requires(T t, double q) {
            {t.evaluate(q)};
    };
}