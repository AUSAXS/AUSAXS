#include <mini/utility/Landscape.h>

SimpleDataset mini::Landscape::as_dataset() const {
    if (evals.empty()) {throw except::bad_order("Landscape::as_dataset: Cannot get evaluated points before a minimization call has been made.");}
    if (evals.front().vals.size() != 1) {throw except::invalid_operation("Landscape::as_dataset: Only 1D landscapes are convertible to datasets.");}

    unsigned int N = evals.size();
    std::vector<double> x(N), y(N);
    bool ordered = true;

    for (unsigned int i = 0; i < N; i++) {
        x[i] = evals[i].vals[0];
        y[i] = evals[i].fval;

        if (i > 0 && x[i] < x[i - 1]) {ordered = false;}
    }
    auto ds = SimpleDataset(x, y, "x", "f(x)");
    if (!ordered) {ds.sort_x();}
    return ds;
}

std::string mini::Landscape::to_string() const {
    std::string s;

    for (const auto& eval : evals) {
        s += eval.to_string() + "\n";
    }
    return s;
}