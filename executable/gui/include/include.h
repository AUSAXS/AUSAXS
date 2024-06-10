#pragma once

#include <elements.hpp>
#include <nfd.hpp>

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
}

struct ColorManager {
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
		for (auto& func : managed_text) {func(get_text_color());}
		for (auto& func : managed_backgrounds) {func(get_color_background());}
	}

	static void manage_input_box(std::shared_ptr<gui::basic_input_box> box) {
		manage_text([box] (gui::color color) {box->set_color(color);});
	}

	static void manage_text(std::function<void(gui::color)>&& func) {
		managed_text.emplace_back(std::move(func));
	}

	static void manage_background(std::function<void(gui::color)>&& func) {
		managed_backgrounds.emplace_back(std::move(func));
	}

	static auto new_background_color() {
		auto bg = std::make_shared<gui::box_element>(get_color_background());
		managed_backgrounds.emplace_back([bg] (gui::color color) {bg->_color = color;});
		return hold(bg);
	}

	inline static bool dark_mode = true;
	inline static std::list<gui::box_element> bgs;
	inline static std::list<std::function<void(gui::color)>> managed_text;
	inline static std::list<std::function<void(gui::color)>> managed_backgrounds;

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
			return shell::Command(utility::remove_all(res.out, "\n\r"));
		}

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
	get_plotter_cmd().append(path).execute();
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