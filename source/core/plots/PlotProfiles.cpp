// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <plots/PlotProfiles.h>

#include <dataset/SimpleDataset.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs::plots;

PlotProfiles::PlotProfiles(observer_ptr<hist::DistanceHistogram> data, const io::File& path) {
	quick_plot(data, path);
}

PlotProfiles::~PlotProfiles() = default;

void PlotProfiles::quick_plot(observer_ptr<hist::DistanceHistogram> data, const io::File& path) {
	PlotHistogram plot;

	if (auto cast = dynamic_cast<hist::ICompositeDistanceHistogram*>(data)) {
		plot.plot(cast->get_profile_aa(), plots::PlotOptions(
			{{"color", style::color::orange}, {"legend", "aa"}, {"normalize", true}, {"xlabel", "q"}, {"ylabel", "I(q)"}, {"logx", true}, {"logy", true}})
		);
		plot.plot(cast->get_profile_aw(), plots::PlotOptions({{"color", style::color::green},  {"legend", "aw"}, {"normalize", true}}));
		plot.plot(cast->get_profile_ww(), plots::PlotOptions({{"color", style::color::blue},   {"legend", "ww"}, {"normalize", true}}));
	} else {
		return;
	}

	if (auto cast = dynamic_cast<const hist::ICompositeDistanceHistogramExv*>(data)) {
		plot.plot(cast->get_profile_ax(), plots::PlotOptions({{"color", style::color::red},    {"normalize", true}, {"legend", "ax"}}));
		plot.plot(cast->get_profile_xx(), plots::PlotOptions({{"color", style::color::purple}, {"normalize", true}, {"legend", "xx"}}));
		plot.plot(cast->get_profile_wx(), plots::PlotOptions({{"color", style::color::cyan},   {"normalize", true}, {"legend", "wx"}}));

		// extra plot: aa / xx
		auto aa = cast->get_profile_aa().get_counts();
		auto xx = cast->get_profile_xx().get_counts();
		std::vector<double> aa_xx;
		for (unsigned int i = 0; i < aa.size(); ++i) {
			aa_xx.push_back(xx[i] / aa[i]);
		}
		PlotHistogram plot2;
		plot2.plot(hist::Histogram(aa_xx, cast->get_profile_aa().get_axis()), plots::PlotOptions(
			{{"xlabel", "q"}, {"ylabel", "I(q)"}, {"logx", true}, {"logy", true}, {"normalize", true}, {"legend", "xx/aa"}})
		);
		plot2.save(io::File(path.append("_xx_aa")));
	}

	plot.save(path);
}