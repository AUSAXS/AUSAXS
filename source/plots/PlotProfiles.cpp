#include <plots/PlotProfiles.h>

#include <dataset/SimpleDataset.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>

using namespace plots;

PlotProfiles::PlotProfiles(const view_ptr<hist::CompositeDistanceHistogram> data, const io::File& path) {
	quick_plot(data, path);
}

PlotProfiles::~PlotProfiles() = default;

void PlotProfiles::quick_plot(const view_ptr<hist::CompositeDistanceHistogram> data, const io::File& path) {
	PlotHistogram plot;
	plot.plot(data->get_profile_aa(), plots::PlotOptions({{"color", style::color::orange}, {"legend", "aa"}, {"normalize", true}, {"xlabel", "q"}, {"ylabel", "I(q)"}, {"logx", true}, {"logy", true}}));
	plot.plot(data->get_profile_aw(), plots::PlotOptions({{"color", style::color::green},  {"legend", "aw"}, {"normalize", true}}));
	plot.plot(data->get_profile_ww(), plots::PlotOptions({{"color", style::color::blue},   {"legend", "ww"}, {"normalize", true}}));

	// if this is a pointer to a CompositeDistanceHistogramFFAvg or CompositeDistanceHistogramFFExplicit, additionally plot the exv
	if (const hist::CompositeDistanceHistogramFFAvg* data_ff_avg = dynamic_cast<const hist::CompositeDistanceHistogramFFAvg*>(data.get())) {
		plot.plot(data_ff_avg->get_profile_ax(), plots::PlotOptions({{"color", style::color::red},    {"normalize", true}, {"legend", "ax"}}));
		plot.plot(data_ff_avg->get_profile_xx(), plots::PlotOptions({{"color", style::color::purple}, {"normalize", true}, {"legend", "xx"}}));
		plot.plot(data_ff_avg->get_profile_wx(), plots::PlotOptions({{"color", style::color::cyan},   {"normalize", true}, {"legend", "wx"}}));

		// extra plot: aa / xx
		auto aa = data->get_profile_aa().get_counts();
		auto xx = data_ff_avg->get_profile_xx().get_counts();
		std::vector<double> aa_xx;
		for (unsigned int i = 0; i < aa.size(); ++i) {
			aa_xx.push_back(xx[i] / aa[i]);
		}
		PlotHistogram plot2;
		plot2.plot(hist::Histogram(aa_xx, data->get_profile_aa().get_axis()), plots::PlotOptions({{"xlabel", "q"}, {"ylabel", "I(q)"}, {"logx", true}, {"logy", true}, {"normalize", true}, {"legend", "xx/aa"}}));
		plot2.save(io::File(path.append("_xx_aa")));
	}

	plot.save(path);
}