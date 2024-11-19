#pragma once

#include <hist/distance_calculator/IHistogramManager.h>
#include <hist/detail/BodyTracker.h>

namespace ausaxs::hist {
	class IPartialHistogramManager : public IHistogramManager, public BodyTracker {
        public:
            IPartialHistogramManager(observer_ptr<const data::Molecule> protein) : BodyTracker(protein) {}
            virtual ~IPartialHistogramManager() override = default;
    };
}