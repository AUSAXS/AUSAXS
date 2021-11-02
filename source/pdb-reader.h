#pragma once

// includes
#include <string>

// my own includes
#include "Reader.h"

class pdb_reader : Reader {
public: 
    std::vector<Atom> read() override {
    };

private:
    Atom read_line() override {
    };
};