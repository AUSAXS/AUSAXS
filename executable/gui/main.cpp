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

auto io_menu(gui::view& view) {
	static auto saxs_box_bg = gui::box(bg_color);
	static auto map_box_bg = gui::box(bg_color);
	static auto output_box_bg = gui::box(bg_color);
	static auto saxs_box = gui::input_box("pdb path");
	static auto map_box = gui::input_box("map path");
	static auto output_box = gui::input_box("output path");
	static bool map_ok = false;
	static bool saxs_ok = false;
	static bool default_output = true;
	output_box.second->set_text("output/em_fitter");

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
			if (constants::filetypes::structure.validate(tmp)) {
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
			if (constants::filetypes::structure.validate(tmp)) {
				matches.push_back(tmp.path());
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
		if (!constants::filetypes::structure.validate(file)) {
			std::cout << "invalid saxs file " << file.path() << std::endl;
			saxs_box_bg = bred;
			saxs_ok = false;
			return;
		}

		settings::saxs_file = file.path();
		std::cout << "saxs file was set to " << settings::saxs_file << std::endl;
		saxs_box_bg = bgreen;
		saxs_ok = true;

		if (map_ok && default_output) {
			std::string path = "output/em_fitter/" + io::File(settings::map_file).stem() + "/" + io::File(settings::saxs_file).stem();
			output_box.second->set_text(path);
			output_box.second->on_enter(path);
		}
	};

	output_box.second->on_enter = [] (std::string_view text) {
		default_output = false;
		settings::general::output = text;
		std::cout << "output path was set to " << settings::general::output << std::endl;
	};

	return gui::htile(
		gui::align_left(
			gui::layer(
				gui::hgrid({0}, map_box.first),
				link(map_box_bg)
			)
		),
		gui::align_center(
			gui::layer(
				gui::hgrid({0}, saxs_box.first),
				link(saxs_box_bg)
			)
		),
		gui::align_right(gui::hgrid({0}, output_box.first))
	);
}

int main(int argc, char* argv[]) {
	gui::app app(argc, argv, "EM fitter", "com.saxs.gui");
	gui::window win(app.name(), 15, gui::rect{20, 20, 1600, 1000});
	win.on_close = [&app]() {app.stop();};

	gui::view view(win);
	auto background = gui::box(bg_color);
	auto logo = gui::scale_element(0.2, gui::image(abs_path("temp/logo.png").c_str()));
	auto settings = gui::vgrid(
		{0},
		io_menu(view)
	);

	view.content
	(
		settings,
		gui::align_right_top(
			gui::margin({10, 10, 10, 10}, logo)
		),
		background
	);

	app.run();
	return 0;
}
