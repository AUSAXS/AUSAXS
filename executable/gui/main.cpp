/*=============================================================================
	Copyright (c) 2016-2020 Joel de Guzman

	Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/
#include <elements.hpp>

#include <data/Molecule.h>
#include <em/ImageStack.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <plots/PlotDistance.h>
#include <plots/PlotIntensity.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotIntensityFitResiduals.h>
#include <constants/Constants.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/All.h>
#include <utility/Limit2D.h>

#include <list>
#include <bitset>
#include <algorithm>
#include <random>
#include <iostream>
#include <memory>
#include <thread>
#include <chrono>

namespace gui = cycfi::elements;

auto constexpr bg_color_accent = gui::rgba(55, 55, 57, 255);
auto constexpr bg_color = gui::rgba(35, 35, 37, 255);
auto constexpr bgreen   = gui::colors::green.level(0.7).opacity(0.4);
auto constexpr bred     = gui::colors::red.level(0.7).opacity(0.4);

auto abs_path(const std::string& path) {
	return std::filesystem::current_path().string() + "/" + path;
}
auto perform_plot(const std::string& path) {
	std::string command = "python3 scripts/plot.py " + path;
	std::system(command.c_str());
};

namespace settings {
	std::string map_file, saxs_file;
	bool lock = false;
}

namespace setup {
	std::unique_ptr<SimpleDataset> saxs_dataset;
	std::unique_ptr<em::ImageStack> map;
}

auto io_menu(gui::view& view) {
	static auto saxs_box_bg = gui::box(bg_color);
	static auto map_box_bg = gui::box(bg_color);
	static auto output_box_bg = gui::box(bg_color);
	static auto saxs_box = gui::input_box("saxs path");
	static auto map_box = gui::input_box("map path");
	static auto output_box = gui::input_box("output path");
	static bool default_output = true;
	output_box.second->set_text("output/em_fitter");
	static bool map_ok = false, saxs_ok = false;

	map_box.second->on_text = [&view] (std::string_view text) {
		static unsigned int last_size = 0;
		if (text.size() == 1) {
			map_box_bg = bg_color_accent;
		} else if (text.empty()) {
			map_box_bg = bg_color;
		}

		if (map_ok) {
			map_ok = false;
			map_box_bg = bg_color_accent;
		}

		// prevent autocompletion when deleting text
		if (text.size() < last_size) {
			last_size = text.size();
			return;
		}
		last_size = text.size();

		// only autocomplete if the last character is a '/' and there are less than 20 matches
		if (text.back() != '/') {return;}
		if (20 < std::distance(std::filesystem::directory_iterator(text), std::filesystem::directory_iterator{})) {return;}

		std::list<std::string> matches;
		for (auto& p : std::filesystem::directory_iterator(text)) {
			io::File tmp(p.path().string());
			if (constants::filetypes::em_map.validate(tmp)) {
				matches.push_back(tmp.path());
			}
		}
		if (matches.empty()) {return;}
		if (matches.size() == 1) {
			// only one match, auto-fill
			settings::map_file = matches.front();
			map_box.second->set_text(matches.front());
			map_box.second->on_enter(matches.front());
			return;
		}
		// find the longest common prefix
		std::string prefix = matches.front();
		for (auto& match : matches) {
			if (prefix == match) {continue;}
			std::string tmp;
			for (size_t i = 0; i < std::min(prefix.size(), match.size()); ++i) {
				if (prefix[i] == match[i]) {
					tmp += prefix[i];
				} else {
					break;
				}
			}
			prefix = tmp;
		}
		if (prefix.size() > 1) {
			map_box.second->set_text(prefix);
		}
	};

	map_box.second->on_enter = [&view] (std::string_view text) {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::em_map.validate(file)) {
			std::cout << "invalid map file " << file.path() << std::endl;
			map_box_bg = bred;
			map_ok = false;
			return;
		}

		settings::map_file = file.path();
		std::cout << "map file was set to " << settings::map_file << std::endl;
		map_box_bg = bgreen;
		map_ok = true;

		if (20 < std::distance(std::filesystem::directory_iterator(file.directory().path()), std::filesystem::directory_iterator{})) {return;}
		for (auto& p : std::filesystem::directory_iterator(file.directory().path())) {
			io::File tmp(p.path().string());
			if (constants::filetypes::saxs_data.validate(tmp)) {
				settings::saxs_file = tmp.path();
				saxs_box.second->set_text(tmp.path());
				saxs_box.second->on_enter(tmp.path());
			} else if (constants::filetypes::setting.validate(tmp)) {
				std::cout << "discovered settings file " << tmp.path() << std::endl;
				settings::read(tmp);
			}
		}

		if (saxs_ok && default_output) {
			std::string path = "output/em_fitter/" + io::File(settings::map_file).stem() + "/" + io::File(settings::saxs_file).stem();
			output_box.second->set_text(path);
			output_box.second->on_enter(path);
		}
	};

	saxs_box.second->on_text = [&view] (std::string_view text) {
		static unsigned int last_size = 0;
		if (text.size() == 1) {
			saxs_box_bg = bg_color_accent;
		} else if (text.empty()) {
			saxs_box_bg = bg_color;
		}

		if (saxs_ok) {
			saxs_ok = false;
			saxs_box_bg = bg_color_accent;
		}

		// prevent autocompletion when deleting text
		if (text.size() < last_size) {
			last_size = text.size();
			return;
		}
		last_size = text.size();

		// only autocomplete if the last character is a '/' and there are less than 20 matches
		if (text.back() != '/') {return;}
		if (20 < std::distance(std::filesystem::directory_iterator(text), std::filesystem::directory_iterator{})) {return;}

		std::list<std::string> matches;
		for (auto& p : std::filesystem::directory_iterator(text)) {
			io::File tmp(p.path().string());
			if (constants::filetypes::saxs_data.validate(tmp)) {
				matches.push_back(tmp.path());
			} else {
				std::cout << "reject possible match: " << tmp.path() << std::endl;
			}
		}
		if (matches.empty()) {return;}
		if (matches.size() == 1) {
			// only one match, auto-fill
			settings::saxs_file = matches.front();
			saxs_box.second->set_text(matches.front());
			saxs_box.second->on_enter(matches.front());
			return;
		}
		// find the longest common prefix
		std::string prefix = matches.front();
		for (auto& match : matches) {
			if (prefix == match) {continue;}
			std::string tmp;
			for (size_t i = 0; i < std::min(prefix.size(), match.size()); ++i) {
				if (prefix[i] == match[i]) {
					tmp += prefix[i];
				} else {
					break;
				}
			}
			prefix = tmp;
		}
		if (prefix.size() > 1) {
			saxs_box.second->set_text(prefix);
		}
	};

	saxs_box.second->on_enter = [&view] (std::string_view text) {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::saxs_data.validate(file)) {
			std::cout << "invalid saxs file " << file.path() << std::endl;
			saxs_box_bg = bred;
			saxs_ok = false;
			setup::saxs_dataset = nullptr;
			return;
		}

		settings::saxs_file = file.path();
		std::cout << "saxs file was set to " << settings::saxs_file << std::endl;
		saxs_box_bg = bgreen;
		saxs_ok = true;
		setup::saxs_dataset = std::make_unique<SimpleDataset>(settings::saxs_file);

		if (map_ok) {
		 	if (default_output) {
				std::string path = "output/em_fitter/" + io::File(settings::map_file).stem() + "/" + io::File(settings::saxs_file).stem();
				output_box.second->set_text(path);
				output_box.second->on_enter(path);
			}
		}
	};

	output_box.second->on_text = [] (std::string_view) {
		default_output = false;
	};

	output_box.second->on_enter = [] (std::string_view text) {
		settings::general::output = text;
		std::cout << "output path was set to " << settings::general::output << std::endl;
	};

	return gui::htile(
		gui::htile(
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					gui::layer(
						map_box.first,
						link(map_box_bg)
					)
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					gui::layer(
						saxs_box.first,
						link(saxs_box_bg)
					)
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					output_box.first
				)
			)
		)
	);
}

auto q_slider(gui::view& view) {
	static auto track = gui::basic_track<5, false>(gui::colors::black);
	static auto thumb = gui::margin(
		{1, 2, 1, 2},
		gui::box(gui::colors::white_smoke)
	);
	static auto qslider = gui::range_slider(
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		gui::slider_labels<9>(
			gui::slider_marks_log<20, 4>(track), 0.8, "1e-4", "1e-3", "1e-2", "1e-1", "1e0"
		),
		{0.1, 0.5}
	);

	static auto qmin_textbox = gui::input_box("q_min");
	static auto qmax_textbox = gui::input_box("q_max");
	static auto qinfo_box = gui::label("test");
	static auto qmin_bg = gui::box(bg_color);
	static auto qmax_bg = gui::box(bg_color);

	auto pretty_printer = [] (float value) {
		std::stringstream ss;
		ss << std::setprecision(2) << std::scientific  << value;
		return ss.str();
	};

	auto axis_transform_inv = [] (float value) {
		return (std::log10(value)-std::log10(constants::axes::q_axis.min))/(std::log10(constants::axes::q_axis.max)-std::log10(constants::axes::q_axis.min));
	};

	auto axis_transform = [] (float x) {
		double logy = x*(std::log10(constants::axes::q_axis.max)-std::log10(constants::axes::q_axis.min)) + std::log10(constants::axes::q_axis.min);
		return std::pow(10, logy);
	};

	qslider.on_change.first = [&view, pretty_printer, axis_transform] (float value) {
		value = axis_transform(value);
		qmin_textbox.second->set_text(pretty_printer(value));
		if (setup::saxs_dataset) {
			unsigned int removed_elements = 0;
			for (unsigned int i = 0; i < setup::saxs_dataset->size(); ++i) {
				auto x = setup::saxs_dataset->x(i);
				removed_elements += !(value < x && x < qslider.value_second());
			}
			if (removed_elements != 0) {
				qinfo_box.set_text("note: ignoring " + std::to_string(removed_elements) + " lines in SAXS file" 
										+ std::string(std::min<int>(4-std::to_string(removed_elements).size(), 0), ' '));
			} else {
				qinfo_box.set_text("");
			}
		} else {
			qinfo_box.set_text("");
		}
		view.refresh(qmax_textbox.first);
		view.refresh(qinfo_box);
	};

	qslider.on_change.second = [&view, pretty_printer, axis_transform] (float value) {
		value = axis_transform(value);
		qmax_textbox.second->set_text(pretty_printer(value));
		if (setup::saxs_dataset) {
			unsigned int removed_elements = 0;
			for (unsigned int i = 0; i < setup::saxs_dataset->size(); ++i) {
				auto x = setup::saxs_dataset->x(i);
				removed_elements += !(qslider.value_first() < x && x < value);
			}
			if (removed_elements != 0) {
				qinfo_box.set_text("note: ignoring " + std::to_string(removed_elements) + " lines in SAXS file" 
										+ std::string(std::min<int>(4-std::to_string(removed_elements).size(), 0), ' '));
			} else {
				qinfo_box.set_text("");
			}
		} else {
			qinfo_box.set_text("");
		}
		view.refresh(qmin_textbox.first);
		view.refresh(qinfo_box);
	};

	qmin_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmin_bg = bg_color;
		} else {
			qmin_bg = bg_color_accent;
		}
	};

	qmin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			qslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			qmin_bg = bg_color;
			view.refresh(qslider);
		} catch (std::exception&) {
			qmin_bg = bred;
		}
	};

	qmax_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmax_bg = bg_color;
		} else {
			qmax_bg = bg_color_accent;
		}
	};

	qmax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			qslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			qmax_bg = bg_color;
			view.refresh(qslider);
		} catch (std::exception&) {
			qmax_bg = bred;
		}
	};

	return gui::vtile(
		gui::margin(
			{50, 0, 50, 0},
			gui::layer(
				link(qslider)
			)
		),
		gui::layer(
			gui::align_left(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							link(qmin_textbox.first),
							link(qmin_bg)
						)
					)
				)
			),
			gui::align_center(
				gui::margin(
					{50, 10, 50, 10},
					link(qinfo_box)
				)
			),
			gui::align_right(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							link(qmax_textbox.first),
							link(qmax_bg)
						)
					)
				)
			)
		)
	);
}

auto alpha_level_slider(gui::view& view) {
	static auto track = gui::basic_track<5, false>(gui::colors::black);
	static auto aslider = gui::range_slider(
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		gui::slider_labels<11>(
			gui::slider_marks<20, 10*5, 10>(track), 0.8, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"
		),
		{0.05, 0.8}
	);

	static auto amin_textbox = gui::input_box("min level");
	static auto amax_textbox = gui::input_box("max level");
	static auto astep_textbox = gui::input_box("steps");
	static auto amin_bg = gui::box(bg_color);
	static auto amax_bg = gui::box(bg_color);
	static auto astep_bg = gui::box(bg_color);

	auto pretty_printer = [] (float value) {
		std::stringstream ss;
		ss << std::setprecision(3) << value;
		return ss.str();
	};

	auto axis_transform = [] (float value) {
		return value*10;
	};

	auto axis_transform_inv = [] (float value) {
		return value/10;
	};

	aslider.on_change.first = [&view, pretty_printer, axis_transform] (float value) {
		amin_textbox.second->set_text(pretty_printer(axis_transform(value)));
		view.refresh(amin_textbox.first);
	};

	aslider.on_change.second = [&view, pretty_printer, axis_transform] (float value) {
		amax_textbox.second->set_text(pretty_printer(axis_transform(value)));
		view.refresh(amax_textbox.first);
	};

	amin_textbox.second->on_text = [&view] (std::string_view text) {
		if (text.empty()) {
			amin_bg = bg_color;
		} else {
			amin_bg = bg_color_accent;
		}
	};

	amin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			aslider.edit_value_first(axis_transform_inv(std::stof(std::string(text))));
			amin_bg = bg_color;
			view.refresh(aslider);
		} catch (std::exception&) {
			amin_bg = bred;
		}
		view.refresh(aslider);
	};

	amax_textbox.second->on_text = [&view] (std::string_view text) {
		if (text.empty()) {
			amax_bg = bg_color;
		} else {
			amax_bg = bg_color_accent;
		}
	};

	amax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			aslider.edit_value_second(axis_transform_inv(std::stof(std::string(text))));
			amax_bg = bg_color;
			view.refresh(aslider);
		} catch (std::exception&) {
			amax_bg = bred;
		}
		view.refresh(aslider);
	};

	return gui::vtile(
		gui::margin(
			{50, 0, 50, 0},
			gui::layer(
				aslider
			)
		),
		gui::layer(
			gui::align_left(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							link(amin_textbox.first),
							link(amin_bg)
						)
					)
				)
			),
			gui::align_center(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							link(astep_textbox.first),
							link(astep_bg)
						)
					)
				)
			),
			gui::align_right(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							link(amax_textbox.first),
							link(amax_bg)
						)
					)
				)
			)
		)
	);
}

auto make_tip(std::string text) {
	return gui::layer(
		gui::margin({20, 8, 20, 8}, gui::basic_text_box(text)), 
		gui::panel{}
	);
}

auto make_misc_settings() {
	auto hydrate = gui::check_box("Hydrate");
	hydrate.value(true);
	hydrate.on_click = [] (bool value) {
		settings::em::hydrate = value;
	};

	auto hydrate_tt = gui::tooltip(
		hydrate,
		make_tip("Hydrate the dummy structure for each fit iteration. This will usually improve the fit substantially.")
	);

	auto fixed_weights = gui::check_box("Fixed weights");
	fixed_weights.value(false);
	fixed_weights.on_click = [] (bool value) {
		settings::em::fixed_weights = value;
	};

	auto fixed_weights_tt = gui::tooltip(
		fixed_weights,
		make_tip("Use fixed weights instead of map densities. This should only be used for maps with a large amount of noise close to the minimum.")
	);

	static auto frequency_bg = gui::box(bg_color);
	auto frequency = gui::input_box("Sample frequency");
	frequency.second->set_text(std::to_string(settings::em::sample_frequency));

	frequency.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			frequency_bg = bg_color;
		} else {
			frequency_bg = bg_color_accent;
		}
	};

	frequency.second->on_enter = [] (std::string_view text) {
		try {
			settings::em::sample_frequency = std::stof(std::string(text));
			frequency_bg = bg_color;
		} catch (std::exception&) {
			std::cout << "invalid sample frequency input" << std::endl;
			frequency_bg = bred;
		}
	};

	static auto frequency_element = gui::hsize(
		200,
		gui::tooltip(
			gui::htile(
				gui::align_right(
					gui::layer(
						gui::hsize(50, frequency.first),
						link(frequency_bg)
					)
				),
				gui::hspace(10),
				gui::align_left(
					gui::label("Sample frequency")
				)
			),
			make_tip("The frequency of sampling the EM grid. Increasing this value will speed up the fit significantly, but will also reduce the accuracy.")
		)
	);

	return gui::htile(
		gui::margin(
			{50, 10, 50, 10},
			hydrate_tt
		),
		gui::margin(
			{50, 10, 50, 10},
			fixed_weights_tt
		),
		gui::margin(
			{50, 10, 50, 10},
			link(frequency_element)
		)
	);
}

auto make_start_button(gui::view& view) {
	static auto start_button = gui::button("start");
	static auto progress_bar = gui::progress_bar(gui::rbox(gui::colors::black), gui::rbox(bgreen));

	auto progress_bar_layout = share(
		gui::margin(
			{10, 100, 10, 100},
			gui::align_center_middle(
				gui::fixed_size(
					{1000, 30},
					link(progress_bar)
				)
			)
		)
	);

	auto start_button_layout = gui::share(
		gui::margin(
			{10, 100, 10, 100},
			gui::align_center_middle(
				gui::hsize(
					200,
					link(start_button)
				)
			)
		)
	);

	static auto content = gui::hold_any(start_button_layout);

	start_button.on_click = [&view, progress_bar_layout] (bool click) {
		if (!setup::saxs_dataset || !setup::map) {
			std::cout << "no saxs data or map file was provided" << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			return;
		}

		std::cout << "##########################" << std::endl;
		std::cout << "### starting EM fitter ###" << std::endl;
		std::cout << "##########################" << std::endl;
		content = progress_bar_layout;

		view.layout(content);
		view.refresh(content);
	};

	return link(content);
}

int main(int argc, char* argv[]) {
	// set maximum qrange, user will be able to change this later
	settings::axes::qmin = 0;
	settings::axes::qmax = 1;

	gui::app app(argc, argv, "EM fitter", "com.saxs.gui");
	gui::window win(app.name(), std::bitset<4>{"1111"}.to_ulong(), gui::rect{20, 20, 1620, 1020});
	win.on_close = [&app]() {app.stop();};

	gui::view view(win);
	auto background = gui::box(bg_color);
	auto header = gui::align_center_top(
		gui::margin(
			{10, 10, 10, 10},
			gui::label("EM fitter")
				.font_size(50)
		)
	);
	auto footer = gui::htile(
		gui::align_left_bottom(
			gui::margin(
				{10, 10, 10, 10},
				gui::label("v0.0.1")
			)
		),
		gui::align_center_bottom(
			gui::margin(
				{10, 10, 10, 10},
				gui::label("Kristian Lytje & Jan Skov Pedersen")
			)
		),
		gui::align_right_bottom(
			gui::margin(
				{10, 10, 10, 10},
				gui::scale_element(0.15, gui::image(abs_path("temp/logo.png").c_str()))
			)
		)
	);
	auto settings = gui::vtile(
		gui::top_margin(
			10,
			gui::label("Input & output")
		),
		io_menu(view),
		gui::hgrid(
			{0.5, 1.0},
			gui::vtile(
				gui::top_margin(
					10,
					gui::label("q-range")
				),
				q_slider(view)
			),
			gui::vtile(
				gui::top_margin(
					10,
					gui::label("alpha levels")
				),
				alpha_level_slider(view)
			)
		),
		make_misc_settings(),
		make_start_button(view)
	);

	view.content (
		gui::vgrid(
			{0.1, 0.9, 1.0},
			header,
			settings,
			footer
		),
		background
	);

	app.run();
	return 0;
}
