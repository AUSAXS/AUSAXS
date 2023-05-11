#pragma once

#include <data/Body.h>

struct BackupBody {
    BackupBody(const Body& body, unsigned int index) : body(body), index(index) {}
    Body body;
    unsigned int index;
};