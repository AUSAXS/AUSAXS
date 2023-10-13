#pragma once

#include <data/Body.h>

struct BackupBody {
    BackupBody(const data::Body& body, unsigned int index) : body(body), index(index) {}
    data::Body body;
    unsigned int index;
};