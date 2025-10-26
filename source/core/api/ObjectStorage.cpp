// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/ObjectStorage.h>

void deallocate(int object_id, int* status) {
    *status = 1;
    ausaxs::api::ObjectStorage::deregister_object(object_id);
    *status = 0;
}