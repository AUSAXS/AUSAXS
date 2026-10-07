// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_data.h>

#include <api/ObjectStorage.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>

#include <string>

using namespace ausaxs;

int data_read(
    const char* filename,
    int* status
) {return execute_with_catch([&]() {
    auto dataset = SimpleDataset(std::string(filename));
    auto data_id = api::ObjectStorage::register_object(std::move(dataset));
    return data_id;
}, status);}

int data_create(
    double* q, double* I, double* Ierr, int n_points,
    int* status
) {return execute_with_catch([&]() {
    if (n_points <= 0) {throw except::invalid_argument("data_create: at least one data point is required.");}
    auto dataset = SimpleDataset(
        std::vector<double>(q, q + n_points),
        std::vector<double>(I, I + n_points),
        std::vector<double>(Ierr, Ierr + n_points)
    );

    // same q-range restriction the file readers apply
    if (settings::axes::clamp_to_qrange) {
        dataset.limit(0, settings::axes::qmin, settings::axes::qmax);
        if (dataset.empty()) {
            throw except::invalid_argument(
                "data_create: no data points inside the q-range "
                "[" + std::to_string(settings::axes::qmin) + ", " + std::to_string(settings::axes::qmax) + "]."
            );
        }
    }
    return api::ObjectStorage::register_object(std::move(dataset));
}, status);}

namespace {
struct _data_get_data_obj {
    explicit _data_get_data_obj(int size) :
        q(size), I(size), Ierr(size)
    {}
    std::vector<double> q, I, Ierr;
};
}
int data_get_data(
    int object_id,
    double** q, double** I, double** Ierr, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto* dataset = api::ObjectStorage::get_object<SimpleDataset>(object_id);
    if (!dataset) {throw except::invalid_argument("Invalid dataset id: \"" + std::to_string(object_id) + "\"");}
    _data_get_data_obj data(dataset->size());
    for (int i = 0; i < dataset->size(); ++i) {
        data.q[i] = dataset->x(i);
        data.I[i] = dataset->y(i);
        data.Ierr[i] = dataset->yerr(i);
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto* ref = api::ObjectStorage::get_object<_data_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *Ierr = ref->Ierr.data();
    *n_points = dataset->size();
    return data_id;
}, status);}