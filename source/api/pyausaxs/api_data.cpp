// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_data.h>
#include <api/ObjectStorage.h>
#include <dataset/SimpleDataset.h>

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

struct _data_get_data_obj {
    explicit _data_get_data_obj(unsigned int size) :
        q(size), I(size), Ierr(size)
    {}
    std::vector<double> q, I, Ierr;
};
int data_get_data(
    int object_id,
    double** q, double** I, double** Ierr, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(object_id);
    if (!dataset) {ErrorMessage::last_error = "Invalid dataset id: \"" + std::to_string(object_id) + "\""; return -1;}
    _data_get_data_obj data(dataset->size());
    for (unsigned int i = 0; i < dataset->size(); ++i) {
        data.q[i] = dataset->x(i);
        data.I[i] = dataset->y(i);
        data.Ierr[i] = dataset->yerr(i);
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_data_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *Ierr = ref->Ierr.data();
    *n_points = static_cast<int>(dataset->size());
    return data_id;
}, status);}