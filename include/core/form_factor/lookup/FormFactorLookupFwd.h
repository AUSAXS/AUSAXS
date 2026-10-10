// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container2D.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorProduct.h>

namespace ausaxs::form_factor::lookup {
    // A table of form factor products, indexed by active slot.
    using table_t = container::Container2D<FormFactorProduct>;
}