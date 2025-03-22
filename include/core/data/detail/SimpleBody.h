#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>

namespace ausaxs::data::detail {
    struct SimpleBody {
        SimpleBody() = default;
        SimpleBody(std::vector<AtomFF>&& atoms) : atoms(std::move(atoms)) {}
        SimpleBody(std::vector<AtomFF>&& atoms, std::vector<Water>&& waters) : atoms(std::move(atoms)), waters(std::move(waters)) {}
        SimpleBody(const std::vector<AtomFF>& atoms) : atoms(atoms) {}
        SimpleBody(const std::vector<AtomFF>& atoms, const std::vector<Water>& waters) : atoms(atoms), waters(waters) {}

        std::vector<AtomFF> atoms;
        std::vector<Water> waters;
    };
}