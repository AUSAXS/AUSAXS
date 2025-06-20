#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/Parameter.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::detail {
    struct Configuration {
        Configuration();
        Configuration(observer_ptr<Rigidbody> rigidbody, double chi2) noexcept;
        ~Configuration();

        std::vector<parameter::Parameter> parameters;
        double chi2;
    };
}