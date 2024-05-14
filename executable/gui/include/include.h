#pragma once

#include <elements.hpp>
#include <nfd.hpp>

#include <em/ImageStack.h>
#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <shell/Command.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <list>

namespace gui = cycfi::elements;
namespace settings {
    inline std::string map_file, pdb_file, saxs_file, output;
}

namespace setup {
	inline std::unique_ptr<data::Molecule> pdb;
	inline std::unique_ptr<SimpleDataset> saxs_dataset;
	inline std::unique_ptr<em::ImageStack> map;
}

#include <fstream>
inline std::ofstream log_file("log.txt");

struct ColorManager {
	static auto new_background_color() {
		managed_background_colors.emplace_back(gui::box_element(dark_mode ? dark_bg : light_bg));
		return link(managed_background_colors.back());
	}

	static auto get_color_background() {
		return dark_mode ? dark_bg : light_bg;
	}

	static auto get_color_accent() {
		return dark_mode ? dark_accent : light_accent;
	}

	static auto get_text_color() {
		return dark_mode ? dark_txt : light_txt;
	}

	static auto get_color_success() {
		return gui::colors::green.level(0.7).opacity(0.4);
	}

	static auto get_color_fail() {
		return gui::colors::red.level(0.7).opacity(0.4);
	}

	static void switch_mode() {
		dark_mode = !dark_mode;
		for (auto& bg : managed_background_colors) {
			if (bg._color == (dark_mode ? light_bg : dark_bg)) {
				bg._color = get_color_background();
			} else if (bg._color == (dark_mode ? light_accent : dark_accent)) {
				bg._color = get_color_background();
			}
		}

		for (auto& label : managed_input_boxes) {
			label->set_color(get_text_color());
		}
	}

	static void manage_input_box(std::shared_ptr<gui::basic_input_box> label) {
		label->set_color(get_text_color());
		managed_input_boxes.push_back(label);
	}

	inline static bool dark_mode = true;
	inline static std::list<gui::box_element> managed_background_colors;
	inline static std::list<std::shared_ptr<gui::basic_input_box>> managed_input_boxes;

	inline static constexpr auto dark_bg = gui::rgba(35, 35, 37, 255);
	inline static constexpr auto light_bg = gui::rgba(255, 255, 255, 255);
	inline static constexpr auto dark_accent = gui::rgba(55, 55, 57, 255);
	inline static constexpr auto light_accent = gui::rgba(200, 200, 200, 255);
	inline static constexpr auto dark_txt = gui::rgba(215, 215, 215, 255);
	inline static constexpr auto light_txt = gui::rgba(65, 65, 65, 255);
};
static auto background = ColorManager::new_background_color();

inline shell::Command get_plotter_cmd() {
	#if defined(_WIN32)
		// first check if plot.exe is available in the path
		auto res = shell::Command("where.exe plot").mute().execute();
		bool plot_exe_available = res.exit_code == 0;
		if (plot_exe_available) {
			log_file << "plot.exe found in path at \"" << res.out << "\"" << std::endl;
			return shell::Command(res.out);
		}
		log_file << "plot.exe not found in path" << std::endl;

		// if not, check if python & the python script is available
		bool python_available = shell::Command("python --version").mute().execute().exit_code == 0;
		bool python_script_available = io::File("scripts/plot.py").exists();
		if (python_available && python_script_available) {
			return shell::Command("python scripts/plot.py");
		}
	#elif defined(__linux__)
		// check if python & the python script is available
		bool python_available = shell::Command("python3 --version").mute().execute().exit_code == 0;
		bool python_script_available = io::File("scripts/plot.py").exists();
		if (python_available && python_script_available) {
			return shell::Command("python3 scripts/plot.py");
		}
	#elif defined (__APPLE__)
		throw std::runtime_error("macOS is not currently supported");
	#endif

	console::print_warning("No plotting utility was found. Please ensure the plot executable is available in the current directory or system path, or that python is installed and the plot.py script is available in the script/ directory.");
	throw std::runtime_error("No plotting utility was found. Please ensure the plot executable is available in the current directory or system path, or that python is installed and the plot.py script is available in the script/ directory.");
}

inline auto perform_plot(const std::string& path) {
	auto cmd = get_plotter_cmd().append(path);
	log_file << "PERFORMING PLOTTER CMD: \"" << cmd.get() << std::endl;
	cmd.execute();
};

enum class NFD_TARGET {FILE, FOLDER};
template<NFD_TARGET target>
inline auto make_dialog_button = [] (auto& text_field, auto& bg, std::pair<std::string, std::string> filter) {
	auto folder_button = gui::button("");
	auto clear_button = gui::button("");
	auto folder_icon = gui::icon(gui::icons::folder_open_empty);
	auto clear_icon = gui::icon(gui::icons::cancel);

	folder_button.on_click = [&text_field, filter] (bool) {
		NFD::Guard guard;
		NFD::UniquePath output;
		nfdresult_t result;

		if constexpr (target == NFD_TARGET::FILE) {
			nfdfilteritem_t filterItem[1] = {{filter.first.c_str(), filter.second.c_str()}};
			result = NFD::OpenDialog(output, filterItem, 1);
		} else {
			result = NFD::PickFolder(output);
		}

		if (result == NFD_OKAY) {
			std::cout << "User picked target: " << output.get() << std::endl;
			text_field.second->set_text(output.get());
			text_field.second->on_enter(output.get());
		} else if (result == NFD_CANCEL) {
			puts("User cancelled selection.");
		} else {
			printf("Error: %s\n", NFD_GetError());
		}
	};

	clear_button.on_click = [&text_field] (bool) {
		text_field.second->set_text("");
		text_field.second->on_text("");
	};

	return gui::htile(
		gui::fixed_size(
			{ 30, 30 },
			gui::layer(
				gui::align_center_middle(
					folder_icon
				),
				folder_button
			)
		),
		gui::hspace(5),
		gui::layer(
			gui::align_right_middle(
				gui::fixed_size(
					{ 18, 18 },
					gui::layer(
						gui::align_center_middle(
							clear_icon
						),
						clear_button
					)
				)
			),
			link(text_field.first),
			link(bg)
		)
	);
};

inline auto make_file_dialog_button(auto& text_field, auto& bg, std::pair<std::string, std::string> filter) {return make_dialog_button<NFD_TARGET::FILE>(text_field, bg, filter);}
inline auto make_folder_dialog_button(auto& text_field, auto& bg) {return make_dialog_button<NFD_TARGET::FOLDER>(text_field, bg, {});}

inline auto make_tip(std::string text) {
	return gui::layer(
		gui::margin({20, 8, 20, 8}, gui::basic_text_box(text)), 
		gui::panel{}
	);
}

inline auto q_slider(gui::view& view) {
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
		{0, 1}
	);

	static auto qmin_textbox = gui::input_box("q_min");
	static auto qmax_textbox = gui::input_box("q_max");
	static auto qinfo_box = gui::label("");
	static auto qmin_bg = ColorManager::new_background_color();
	static auto qmax_bg = ColorManager::new_background_color();
	ColorManager::manage_input_box(qmin_textbox.second);
	ColorManager::manage_input_box(qmax_textbox.second);

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

	auto update_removed_counter = [&view] () {
		static bool extend_wait = false;
		static bool waiting = false;
		static std::thread worker;

		// extend wait period if the user is still moving the slider
		if (waiting) {
			extend_wait = true;
			return;	
		}

		// a valid dataset must be present
		if (!setup::saxs_dataset) {return;}

		// ensure that the worker thread is finished and can be destructed
		if (worker.joinable()) {
			worker.join();
		}

		waiting = true;
		worker = std::thread([&view] () {
			do {
				extend_wait = false;
				std::this_thread::sleep_for(std::chrono::milliseconds(500));			
			} while (extend_wait);

			unsigned int removed_elements = 0;
			for (unsigned int i = 0; i < setup::saxs_dataset->size(); ++i) {
				auto x = setup::saxs_dataset->x(i);
				removed_elements += !(qslider.value_first() < x && x < qslider.value_second());
			}
			if (removed_elements != 0) {
				qinfo_box.set_text("note: ignoring " + std::to_string(removed_elements) + " lines in SAXS file" 
										+ std::string(std::min<int>(4-std::to_string(removed_elements).size(), 0), ' '));
			} else {
				qinfo_box.set_text("");
			}
			view.refresh(qinfo_box);
			waiting = false;
		});
	};

	qslider.on_change.first = [&view, pretty_printer, axis_transform, update_removed_counter] (float value) {
		value = axis_transform(value);
		qmin_textbox.second->set_text(pretty_printer(value));
		update_removed_counter();
		view.refresh(qmin_textbox.first);
		view.refresh(qinfo_box);
	};

	qslider.on_change.second = [&view, pretty_printer, axis_transform, update_removed_counter] (float value) {
		value = axis_transform(value);
		qmax_textbox.second->set_text(pretty_printer(value));
		update_removed_counter();
		view.refresh(qmax_textbox.first);
		view.refresh(qinfo_box);
	};

	qmin_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmin_bg.get() = ColorManager::get_color_background();
		} else {
			qmin_bg.get() = ColorManager::get_color_accent();
		}
	};

	qmin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			qslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			qmin_bg.get() = ColorManager::get_color_background();
			view.refresh(qslider);
		} catch (std::exception&) {
			qmin_bg.get() = ColorManager::get_color_fail();
		}
		return true;
	};

	qmax_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmax_bg.get() = ColorManager::get_color_background();
		} else {
			qmax_bg.get() = ColorManager::get_color_accent();
		}
	};

	qmax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			qslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			qmax_bg.get() = ColorManager::get_color_background();
			view.refresh(qslider);
		} catch (std::exception&) {
			qmax_bg.get() = ColorManager::get_color_fail();
		}
		return true;
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

/**
 * @brief Text autocompleter. 
 * 
 * First determines if the path is valid for autocompletion.
 * If so, it evaluates all files in the directory @a path with the @a cmp_func. 
 * If there are no matches, an empty string is returned along with 'false'. 
 * If there is only a single match, the full path to the file is returned along with 'true'. 
 * If there are multiple matches, the longest common prefix is returned along with 'false'. 
 */
inline auto autocomplete = [] (std::string_view path, unsigned int& last_size, std::function<bool(const io::File&)> cmp_func) {
	// prevent autocompletion when deleting text
	if (path.size() < last_size) {
		last_size = path.size();
		return std::make_pair(std::string(), false);
	}
	last_size = path.size();

	// only autocomplete if the last character is a '/' and there are less than 20 matches
	if (path.back() != '/') {return std::make_pair(std::string(), false);;}
	if (!std::filesystem::is_directory(path)) {return std::make_pair(std::string(), false);;}
	if (20 < std::distance(std::filesystem::directory_iterator(path), std::filesystem::directory_iterator{})) {return std::make_pair(std::string(), false);;}

	std::list<std::string> matches;
	for (auto& p : std::filesystem::directory_iterator(path)) {
		io::File tmp(p.path().string());
		if (cmp_func(tmp)) {
			matches.push_back(tmp.path());
		}
	}

	// no matches, return empty string
	if (matches.empty()) {return std::make_pair(std::string(path), false);}

	// only one match, auto-fill
	if (matches.size() == 1) {
		return std::make_pair(matches.front(), true);
	}

	// multiple matches, find the longest common prefix
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

	// no slashes
	if (prefix.size() > 1) {
		return std::make_pair(prefix, false);
	}
	return std::make_pair(std::string(path), false);
};

inline auto alpha_level_slider(gui::view& view) {
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
			gui::slider_marks_lin<20, 10, 5>(track), 0.8, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"
		),
		{0.05, 0.8}
	);

	static auto amin_textbox = gui::input_box("min level");
	static auto amax_textbox = gui::input_box("max level");
	static auto amin_infobox = gui::label("          ");
	static auto amax_infobox = gui::label("          ");
	static auto astep_textbox = gui::input_box("steps");
	static auto amin_bg = ColorManager::new_background_color();
	static auto amax_bg = ColorManager::new_background_color();
	static auto astep_bg = ColorManager::new_background_color();

	auto astep_tt = gui::tooltip(
		link(astep_textbox.first),
		make_tip("The requested number of fit evaluations. This is only used as a guideline.")
	);

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

	auto update_mass_range = [&view, axis_transform] () {
		static std::pair<double, double> values = {-1, -1};
		static bool extend_wait = false;
		static bool waiting = false;
		static std::thread worker;

		// extend wait period if the user is still moving the slider
		if (waiting) {
			extend_wait = true;
			return;	
		}

		// a valid map must be present
		if (!setup::map) {return;}

		// ensure that the worker thread is finished and can be destructed
		if (worker.joinable()) {
			worker.join();
		}

		waiting = true;
		worker = std::thread([&view, axis_transform] () {
			do {
				extend_wait = false;
				std::this_thread::sleep_for(std::chrono::milliseconds(500));			
			} while (extend_wait);

			if (aslider.value_first() != values.first) {
				values.first = aslider.value_first();
				double mass_min = setup::map->get_mass(setup::map->from_level(axis_transform(values.first)));
				std::stringstream ss;
				ss << std::setprecision(3) << mass_min;
				amin_infobox.set_text("mass: " + ss.str() + " kDa");
				view.refresh(amin_infobox);
				view.refresh(amin_textbox.first);
			}

			if (aslider.value_second() != values.second) {
				values.second = aslider.value_second();
				double mass_max = setup::map->get_mass(setup::map->from_level(axis_transform(values.second)));
				std::stringstream ss;
				ss << std::setprecision(3) << mass_max;
				amax_infobox.set_text("mass: " + ss.str() + " kDa");
				view.refresh(amax_infobox);
				view.refresh(amax_textbox.first);
			}

			waiting = false;
			return;
		});
	};

	aslider.on_change.first = [&view, pretty_printer, axis_transform, update_mass_range] (float value) {
		amin_textbox.second->set_text(pretty_printer(axis_transform(value)));
		update_mass_range();
		view.refresh(amin_textbox.first);
	};

	aslider.on_change.second = [&view, pretty_printer, axis_transform, update_mass_range] (float value) {
		amax_textbox.second->set_text(pretty_printer(axis_transform(value)));
		update_mass_range();
		view.refresh(amax_textbox.first);
	};

	amin_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {			
			amin_bg.get() = ColorManager::get_color_background();
		} else {
			amin_bg.get() = ColorManager::get_color_accent();
		}
	};

	amin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			aslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			amin_bg.get() = ColorManager::get_color_background();
			view.refresh(aslider);
		} catch (std::exception&) {
			amin_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(aslider);
		return true;
	};

	amax_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			amax_bg.get() = ColorManager::get_color_background();
		} else {
			amax_bg.get() = ColorManager::get_color_accent();
		}
	};

	amax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			aslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			amax_bg.get() = ColorManager::get_color_background();
			view.refresh(aslider);
		} catch (std::exception&) {
			amax_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(aslider);
		return true;
	};

	astep_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			astep_bg.get() = ColorManager::get_color_background();
		} else {
			astep_bg.get() = ColorManager::get_color_accent();
		}
	};

	astep_textbox.second->on_enter = [&view] (std::string_view text) -> bool {
		try {
			settings::fit::max_iterations = std::stof(std::string(text));
			astep_bg.get() = ColorManager::get_color_background();
		} catch (std::exception&) {
			astep_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(astep_textbox.first);
		return true;
	};

	return gui::vtile(
		gui::margin(
			{50, 0, 50, 0},
			gui::layer(
				link(aslider)
			)
		),
		gui::layer(
			gui::align_left(
				gui::margin(
					{50, 10, 50, 10},
					gui::vtile(
						gui::hsize(
							100,
							gui::layer(
								link(amin_textbox.first),
								link(amin_bg)
							)
						),
						gui::margin_top(
							2,
							link(amin_infobox)
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
							astep_tt,
							link(astep_bg)
						)
					)
				)
			),
			gui::align_right(
				gui::margin(
					{50, 10, 50, 10},
					gui::vtile(
						gui::hsize(
							100,
							gui::layer(
								link(amax_textbox.first),
								link(amax_bg)
							)
						),
						gui::margin_top(
							2,
							link(amax_infobox)
						)
					)
				)
			)
		)
	);
}
