#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/detail/BodyTracker.h>

namespace ausaxs::hist {
    template<bool use_weighted_distribution> 
	class IPartialHistogramManager : public HistogramManager<use_weighted_distribution>, public BodyTracker {
        public:
            IPartialHistogramManager(observer_ptr<const data::Molecule> protein) : HistogramManager<use_weighted_distribution>(protein), BodyTracker(protein) {}
            virtual ~IPartialHistogramManager() override = default;
    };
}