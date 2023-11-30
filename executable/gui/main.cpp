/*=============================================================================
	Copyright (c) 2016-2020 Joel de Guzman

	Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/
#include <elements.hpp>
#include <algorithm>
#include <random>
#include <iostream>
#include <memory>

#include <data/Molecule.h>
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
}

namespace setup {
	SimpleDataset saxs_dataset;
}

auto io_menu(gui::view& view) {
	static auto saxs_box_bg = gui::box(bg_color);
	static auto map_box_bg = gui::box(bg_color);
	static auto output_box_bg = gui::box(bg_color);
	static auto saxs_box = gui::input_box("pdb path");
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
			setup::saxs_dataset = SimpleDataset();
			return;
		}

		settings::saxs_file = file.path();
		std::cout << "saxs file was set to " << settings::saxs_file << std::endl;
		saxs_box_bg = bgreen;
		saxs_ok = true;
		setup::saxs_dataset = SimpleDataset(settings::saxs_file);

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
						gui::hgrid({0}, map_box.first),
						link(map_box_bg)
					)
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					gui::layer(
						gui::hgrid({0}, saxs_box.first),
						link(saxs_box_bg)
					)
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					gui::hgrid({0}, output_box.first)
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
	static auto qmin_slider = gui::slider(
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		gui::slider_labels<9>(
			gui::slider_marks<20, 8*5, 8>(track), 0.8, "1e-4", "5e-4", "1e-3", "5e-3", "1e-2", "5e-2", "1e-1", "5e-1", "1e0"
		),
		0.1
	);

	static auto qmax_slider = gui::slider(
		gui::fixed_size(
			{5, 30},
			gui::box(gui::colors::light_gray)
		),
		track,
		0.9
	);

	static auto qmin_textbox = gui::input_box("q_min");
	static auto qmax_textbox = gui::input_box("q_max");
	static auto qinfo_box = gui::label("test");
	static auto qmin_bg = gui::box(bg_color);
	static auto qmax_bg = gui::box(bg_color);

	qmin_slider.on_change = [&view] (float value) {
		qmin_textbox.second->set_text(std::to_string(value));
		if (!setup::saxs_dataset.empty()) {
			unsigned int removed_elements = 0;
			for (; removed_elements < setup::saxs_dataset.size(); ++removed_elements) {
				if (value < setup::saxs_dataset.x(removed_elements) || setup::saxs_dataset.x(removed_elements) < qmax_slider.value()) {
					break;
				}
			}
			qinfo_box.set_text("note: skipping " + std::to_string(removed_elements) + " elements");
		}
		view.refresh(qmin_textbox.first);
		view.refresh(qinfo_box);
	};

	qmax_slider.on_change = [&view] (float value) {
		qmax_textbox.second->set_text(std::to_string(value));
		if (!setup::saxs_dataset.empty()) {
			unsigned int removed_elements = 0;
			for (; removed_elements < setup::saxs_dataset.size(); ++removed_elements) {
				if (qmin_slider.value() < setup::saxs_dataset.x(removed_elements) || setup::saxs_dataset.x(removed_elements) < value) {
					break;
				}
			}
			qinfo_box.set_text("note: skipping " + std::to_string(removed_elements) + " elements");
		}
		view.refresh(qmax_textbox.first);
		view.refresh(qinfo_box);
	};

	qmin_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmin_bg = bg_color;
		} else {
			qmin_bg = bg_color_accent;
		}
	};

	qmin_textbox.second->on_enter = [&view] (std::string_view text) {
		try {
			qmin_slider.edit_value(std::stof(std::string(text)));
			qmin_bg = bg_color;
			view.refresh(qmin_slider);
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

	qmax_textbox.second->on_enter = [&view] (std::string_view text) {
		try {
			qmax_slider.edit_value(std::stof(std::string(text)));
			qmax_bg = bg_color;
			view.refresh(qmax_slider);
		} catch (std::exception&) {
			qmax_bg = bred;
		}
	};

	return gui::vtile(
		gui::margin(
			{50, 0, 50, 0},
			gui::layer(
				qmin_slider,
				qmax_slider
			)
		),
		gui::layer(
			gui::align_left(
				gui::margin(
					{50, 10, 50, 10},
					gui::hsize(
						100,
						gui::layer(
							gui::hgrid({0}, link(qmin_textbox.first)),
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
							gui::hgrid({0}, link(qmax_textbox.first)),
							link(qmax_bg)
						)
					)
				)
			)
		)
	);
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
				gui::label("Kristian Lytje & Jan Petersen")
			)
		),
		gui::align_right_bottom(
			gui::margin(
				{10, 10, 10, 10},
				gui::scale_element(0.15, gui::image(abs_path("temp/logo.png").c_str()))
			)
		)
	);
	auto settings = gui::vgrid(
		{0},
		gui::vtile(
			gui::top_margin(
				10,
				gui::label("Input & output")
			),
			io_menu(view),
			gui::top_margin(
				10,
				gui::label("q-range")
			),
			q_slider(view)
		)
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
