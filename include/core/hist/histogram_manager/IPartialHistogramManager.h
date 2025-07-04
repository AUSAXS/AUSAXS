// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/detail/BodyTracker.h>

namespace ausaxs::hist {
	class IPartialHistogramManager : public IHistogramManager, public BodyTracker {
        public:
            IPartialHistogramManager(observer_ptr<const data::Molecule> protein) : BodyTracker(protein) {}
            virtual ~IPartialHistogramManager() override = default;
    };
}