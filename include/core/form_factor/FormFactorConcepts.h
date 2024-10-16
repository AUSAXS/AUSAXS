#include <concepts>

template<typename T>
    concept FormFactorType = requires(T t, double q) {
        {t.evaluate(q)};
};