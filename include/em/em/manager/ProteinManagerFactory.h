// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/manager/ProteinManager.h>

#include <memory>

namespace ausaxs::em::factory {
    std::unique_ptr<em::managers::ProteinManager> create_manager(observer_ptr<const ImageStackBase> images);
}