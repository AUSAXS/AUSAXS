#pragma once

#include <elements.hpp>
#include <nfd.hpp>

#include <em/ImageStack.h>
#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <shell/Command.h>
#include <constants/Constants.h>
#include <settings/All.h>

#include <list>

namespace gui = cycfi::elements;
namespace settings {
    std::string map_file, pdb_file, saxs_file, output;
}

namespace setup {
	std::unique_ptr<data::Molecule> pdb;
	std::unique_ptr<SimpleDataset> saxs_dataset;
	std::unique_ptr<em::ImageStack> map;
}

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
		return gui::rgba(117, 117, 117, 255);
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
	}

	inline static bool dark_mode = true;
	inline static std::list<gui::box_element> managed_background_colors;

	inline static constexpr auto dark_bg = gui::rgba(35, 35, 37, 255);
	inline static constexpr auto light_bg = gui::rgba(255, 255, 255, 255);
	inline static constexpr auto dark_accent = gui::rgba(55, 55, 57, 255);
	inline static constexpr auto light_accent = gui::rgba(200, 200, 200, 255);
};
static auto background = ColorManager::new_background_color();

shell::Command get_plotter_cmd() {
	#ifdef _WIN32
		// first check if plot.exe is available in the path
		auto res = shell::Command("where plot").mute().execute();
		bool plot_exe_available = res.exit_code == 0;
		if (plot_exe_available) {
			return shell::Command(res.out);
		}

		// if not, check if python & the python script is available
		bool python_available = shell::Command("python --version").mute().execute();
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

	throw std::runtime_error("No plotting utility was found. Please ensure the plot executable is available in the current directory or system path, or that python is installed and the plot.py script is available in the script/ directory.");
}

auto perform_plot(const std::string& path) {
	get_plotter_cmd().append(path).execute();
};

enum class NFD_TARGET {FILE, FOLDER};
template<NFD_TARGET target>
auto make_dialog_button = [] (auto& text_field, auto& bg, std::pair<std::string, std::string> filter) {
   auto button = gui::button("");
   auto icon = gui::icon(gui::icons::folder_open_empty);

   button.on_click = [&text_field, filter] (bool) {
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

   return gui::htile(
      gui::fixed_size(
         { 30, 30 },
         gui::layer(
            gui::align_center_middle(
               icon
            ),
            button
         )
      ),
      gui::hspace(5),
	  gui::layer(
	      link(text_field.first),
		  link(bg)
	  )
   );
};

auto make_file_dialog_button(auto& text_field, auto& bg, std::pair<std::string, std::string> filter) {return make_dialog_button<NFD_TARGET::FILE>(text_field, bg, filter);}
auto make_folder_dialog_button(auto& text_field, auto& bg) {return make_dialog_button<NFD_TARGET::FOLDER>(text_field, bg, {});}

auto make_tip(std::string text) {
	return gui::layer(
		gui::margin({20, 8, 20, 8}, gui::basic_text_box(text)), 
		gui::panel{}
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
		{0, 1}
	);

	static auto qmin_textbox = gui::input_box("q_min");
	static auto qmax_textbox = gui::input_box("q_max");
	static auto qinfo_box = gui::label("");
	static auto qmin_bg = ColorManager::new_background_color();
	static auto qmax_bg = ColorManager::new_background_color();

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

	qmin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			qslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			qmin_bg.get() = ColorManager::get_color_background();
			view.refresh(qslider);
		} catch (std::exception&) {
			qmin_bg.get() = ColorManager::get_color_fail();
		}
	};

	qmax_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			qmax_bg.get() = ColorManager::get_color_background();
		} else {
			qmax_bg.get() = ColorManager::get_color_accent();
		}
	};

	qmax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			qslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			qmax_bg.get() = ColorManager::get_color_background();
			view.refresh(qslider);
		} catch (std::exception&) {
			qmax_bg.get() = ColorManager::get_color_fail();
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

	amin_textbox.second->on_text = [&view] (std::string_view text) {
		if (text.empty()) {			
			amin_bg.get() = ColorManager::get_color_background();
		} else {
			amin_bg.get() = ColorManager::get_color_accent();
		}
	};

	amin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			aslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			amin_bg.get() = ColorManager::get_color_background();
			view.refresh(aslider);
		} catch (std::exception&) {
			amin_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(aslider);
	};

	amax_textbox.second->on_text = [&view] (std::string_view text) {
		if (text.empty()) {
			amax_bg.get() = ColorManager::get_color_background();
		} else {
			amax_bg.get() = ColorManager::get_color_accent();
		}
	};

	amax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) {
		try {
			aslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			amax_bg.get() = ColorManager::get_color_background();
			view.refresh(aslider);
		} catch (std::exception&) {
			amax_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(aslider);
	};

	astep_textbox.second->on_text = [&view] (std::string_view text) {
		if (text.empty()) {
			astep_bg.get() = ColorManager::get_color_background();
		} else {
			astep_bg.get() = ColorManager::get_color_accent();
		}
	};

	astep_textbox.second->on_enter = [&view] (std::string_view text) {
		try {
			settings::fit::max_iterations = std::stof(std::string(text));
			astep_bg.get() = ColorManager::get_color_background();
		} catch (std::exception&) {
			astep_bg.get() = ColorManager::get_color_fail();
		}
		view.refresh(astep_textbox.first);
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
						gui::top_margin(
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
						gui::top_margin(
							2,
							link(amax_infobox)
						)
					)
				)
			)
		)
	);
}
