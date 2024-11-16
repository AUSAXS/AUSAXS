#pragma once

#include <data/Body.h>

namespace ausaxs::rigidbody::transform {
    struct BackupBody {
        BackupBody(const data::Body& body, unsigned int index) : body(body), index(index) {}
        data::Body body;
        unsigned int index;
    };
}