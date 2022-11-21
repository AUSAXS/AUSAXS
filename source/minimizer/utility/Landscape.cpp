#include <mini/utility/Landscape.h>

using namespace mini;

SimpleDataset Landscape::as_dataset() const {
    std::cout << "\tcheckpoint" << std::endl;
    if (evals.empty()) {throw except::bad_order("Landscape::as_dataset: Cannot get evaluated points before a minimization call has been made.");}
    if (evals[0].vals.size() != 1) {throw except::invalid_operation("Landscape::as_dataset: Only 1D landscapes are convertible to datasets.");}

    std::cout << "\tcheckpoint" << std::endl;
    unsigned int N = evals.size();
    std::vector<double> x(N), y(N);
    bool ordered = true;
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evals[i].vals[0];
        y[i] = evals[i].fval;

        if (i > 0 && x[i] < x[i - 1]) {ordered = false;}
    }
    std::cout << "\tcheckpoint" << std::endl;
    auto ds = SimpleDataset(x, y, "x", "f(x)");
    std::cout << "\tcheckpoint" << std::endl;
    if (!ordered) {ds.sort_x();}
    std::cout << "\tcheckpoint" << std::endl;
    return ds;
}