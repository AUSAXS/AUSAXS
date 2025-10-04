// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/Landscape.h>
#include <mini/detail/Evaluation.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

mini::Landscape::Landscape(std::vector<Evaluation>&& evals) : evals(std::move(evals)) {}

mini::Landscape::Landscape(const std::vector<Evaluation>& evals) : evals(evals) {}

void mini::Landscape::append(const std::vector<Evaluation>& evals) {this->evals.insert(this->evals.end(), evals.begin(), evals.end());}

void mini::Landscape::append(const Landscape& evals) {append(evals.evals);}

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