#include <elements.hpp>
#include <nfd.hpp>

#include <filesystem>
#ifdef __APPLE__
static struct _dummy{
	_dummy() {
		std::string out = std::string(std::getenv("HOME")) + "/Documents/ausaxs";
		std::filesystem::create_directories(out);
		std::filesystem::current_path(out);
	}
} _dummy_instance;
#endif

#include <data/Molecule.h>
#include <em/ImageStack.h>
#include <fitter/SmartFitter.h>
#include <plots/PlotDistance.h>
#include <constants/Constants.h>
#include <constants/Version.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <utility/Limit2D.h>
#include <utility/MultiThreading.h>
#include <fitter/FitReporter.h>
#include <shell/Command.h>
#include <constants/ValidFileExtensions.h>
#include <settings/All.h>

#include <gui/helper.h>
#include <gui/logo.h>
#include <gui/resources.h>

#include <iostream>
#include <memory>
#include <thread>
#include <bitset>

namespace gui = cycfi::elements;

namespace setup {
	inline std::unique_ptr<em::ImageStack> map;
}

auto constexpr bg_color_accent = gui::rgba(55, 55, 57, 255);
auto constexpr bg_color = gui::rgba(35, 35, 37, 255);
auto constexpr bgreen   = gui::colors::green.level(0.7).opacity(0.4);
auto constexpr bred     = gui::colors::red.level(0.7).opacity(0.4);
auto plot_names = std::vector<std::pair<std::string, std::string>>{
	{"chi2_evaluated_points_full", "χ²"},
	{"chi2_evaluated_points_limited", "χ² reduced axes"},
	{"chi2_near_minimum", "χ² near minimum"},
	{"chi2_evaluated_points_limited_mass", "χ² reduced axes (mass)"},
	{"chi2_near_minimum_mass", "χ² near minimum (mass)"},
	{"log", "log plot"},
	{"loglog", "log-log plot"}
};

auto abs_path(const std::string& path) {
	return std::filesystem::current_path().string() + "/" + path;
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

	map_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			map_box_bg = bg_color_accent;
		} else if (text.empty()) {
			map_box_bg = bg_color;
		}

		if (map_ok) {
			map_ok = false;
			map_box_bg = bg_color_accent;
		}

		static unsigned int last_size = 0;
		auto fill = autocomplete(text, last_size, [] (const io::File& p) {return constants::filetypes::em_map.check(p);});
		if (!fill.first.empty()) {map_box.second->set_text(fill.first);}
		if (fill.second) {map_box.second->on_enter(fill.first);}
	};

	map_box.second->on_enter = [] (std::string_view text) -> bool {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::em_map.check(file)) {
			std::cout << "invalid map file " << file.path() << std::endl;
			map_box_bg = bred;
			map_ok = false;
			return true;
		}

		settings::map_file = file.path();
		std::cout << "map file was set to " << settings::map_file << std::endl;
		try {
			setup::map = std::make_unique<em::ImageStack>(settings::map_file);
			map_box_bg = bgreen;
			map_ok = true;
		} catch (std::exception& e) {
			map_box_bg = bred;
			map_ok = false;
			setup::map = nullptr;
			return true;
		}

		if (!saxs_ok) {
			if (20 < std::distance(std::filesystem::directory_iterator(file.directory().path()), std::filesystem::directory_iterator{})) {return true;}
			for (auto& p : std::filesystem::directory_iterator(file.directory().path())) {
				io::File tmp(p.path().string());
				if (constants::filetypes::saxs_data.check(tmp)) {
					settings::saxs_file = tmp.path();
					saxs_box.second->set_text(tmp.path());
					saxs_box.second->on_enter(tmp.path());
					break;
				}
			}
		}

		if (saxs_ok && default_output) {
			std::string path = "output/em_fitter/" + io::File(settings::map_file).stem() + "/" + io::File(settings::saxs_file).stem();
			output_box.second->set_text(path);
			output_box.second->on_enter(path);
		}
		return true;
	};

	saxs_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			saxs_box_bg = bg_color_accent;
		} else if (text.empty()) {
			saxs_box_bg = bg_color;
		}

		if (saxs_ok) {
			saxs_ok = false;
			saxs_box_bg = bg_color_accent;
		}

		static unsigned int last_size = 0;
		auto fill = autocomplete(text, last_size, [] (const io::File& p) {return constants::filetypes::saxs_data.check(p);});
		if (!fill.first.empty()) {saxs_box.second->set_text(fill.first);}
		if (fill.second) {saxs_box.second->on_enter(fill.first);}
	};

	saxs_box.second->on_enter = [] (std::string_view text) -> bool {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::saxs_data.check(file)) {
			std::cout << "invalid saxs file " << file.path() << std::endl;
			saxs_box_bg = bred;
			saxs_ok = false;
			setup::saxs_dataset = nullptr;
			return true;
		}

		settings::saxs_file = file.path();
		std::cout << "saxs file was set to " << settings::saxs_file << std::endl;
		try {
			setup::saxs_dataset = std::make_unique<SimpleDataset>(settings::saxs_file);
			saxs_box_bg = bgreen;
			saxs_ok = true;
		} catch (std::exception& e) {
			saxs_box_bg = bred;
			saxs_ok = false;
			setup::saxs_dataset = nullptr;
			return true;
		}

		if (map_ok) {
		 	if (default_output) {
				std::string path = "output/em_fitter/" + io::File(settings::map_file).stem() + "/" + io::File(settings::saxs_file).stem();
				output_box.second->set_text(path);
				output_box.second->on_enter(path);
			}
		}
		return true;
	};

	output_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			output_box_bg = bg_color_accent;
		} else if (text.empty()) {
			output_box_bg = bg_color;
		}
		default_output = false;
	};

	output_box.second->on_enter = [&view] (std::string_view text) -> bool {
		settings::general::output = text;
		if (settings::general::output.back() != '/') {
			settings::general::output += "/";
			view.refresh(output_box.first);
		}
		std::cout << "output path was set to " << settings::general::output << std::endl;
		return true;
	};

	auto map_box_field = make_file_dialog_button(map_box, map_box_bg, {"EM map", "map,ccp4,mrc"});
	auto saxs_box_field = make_file_dialog_button(saxs_box, saxs_box_bg, {"SAXS data", "dat,scat"});
	auto output_box_field = make_folder_dialog_button(output_box, output_box_bg);

	return gui::htile(
		gui::htile(
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					map_box_field
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					saxs_box_field
				)
			),
			gui::margin(
				{50, 10, 50, 10},
				gui::hsize(
					300,
					output_box_field
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
	static auto amin_bg = gui::box(bg_color);
	static auto amax_bg = gui::box(bg_color);
	static auto astep_bg = gui::box(bg_color);

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
			amin_bg = bg_color;
		} else {
			amin_bg = bg_color_accent;
		}
	};

	amin_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			aslider.value_first(axis_transform_inv(std::stof(std::string(text))));
			amin_bg = bg_color;
			view.refresh(aslider);
		} catch (std::exception&) {
			amin_bg = bred;
		}
		view.refresh(aslider);
		return true;
	};

	amax_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			amax_bg = bg_color;
		} else {
			amax_bg = bg_color_accent;
		}
	};

	amax_textbox.second->on_enter = [&view, axis_transform_inv] (std::string_view text) -> bool {
		try {
			aslider.value_second(axis_transform_inv(std::stof(std::string(text))));
			amax_bg = bg_color;
			view.refresh(aslider);
		} catch (std::exception&) {
			amax_bg = bred;
		}
		view.refresh(aslider);
		return true;
	};

	astep_textbox.second->on_text = [] (std::string_view text) {
		if (text.empty()) {
			astep_bg = bg_color;
		} else {
			astep_bg = bg_color_accent;
		}
	};

	astep_textbox.second->on_enter = [&view] (std::string_view text) -> bool {
		try {
			settings::fit::max_iterations = std::stof(std::string(text));
			astep_bg = bg_color;
		} catch (std::exception&) {
			astep_bg = bred;
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

	frequency.second->on_enter = [] (std::string_view text) -> bool {
		try {
			settings::em::sample_frequency = std::stof(std::string(text));
			frequency_bg = bg_color;
		} catch (std::exception&) {
			std::cout << "invalid sample frequency input" << std::endl;
			frequency_bg = bred;
		}
		return true;
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

	auto progress_bar_layout = gui::margin(
		{10, 100, 10, 100},
		gui::align_center_middle(
			gui::fixed_size(
				{1000, 30},
				link(progress_bar)
			)
		)
	);

	auto start_button_layout = gui::margin(
		{10, 100, 10, 100},
		gui::align_center_middle(
			gui::hsize(
				200,
				link(start_button)
			)
		)
	);

	static auto deck = gui::deck_composite();
	deck.push_back(gui::share(start_button_layout));
	deck.push_back(gui::share(progress_bar_layout));

	static std::thread worker;
	start_button.on_click = [&view] (bool) {
		if (!setup::saxs_dataset || !setup::map) {
			std::cout << "no saxs data or map file was provided" << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			return;
		}

		static auto observer = setup::map->get_progress_observer();
		observer->on_notify = [&view] (int progress) {
			progress_bar.value(float(progress)/(2*settings::fit::max_iterations));
			view.refresh(deck);
		};

		deck.select(1);
		view.refresh();
		worker = std::thread([&view] () {
			auto res = setup::map->fit(settings::saxs_file);

			// small animation to make the bar reach 100%
			while (progress_bar.value() < 1) {
				progress_bar.value(progress_bar.value() + 0.02);
				view.refresh(deck);
				std::this_thread::sleep_for(std::chrono::milliseconds(25));
			}

			// perform the plots
			res->curves.save(settings::general::output + "ausaxs.fit", "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof));
			fitter::FitReporter::save(res.get(), settings::general::output + "report.txt");
			perform_plot(settings::general::output);

			auto make_image_pane = [] (const io::File& path) {
				return gui::image(std::filesystem::current_path().string() + "/" + path.path().c_str(), 0.15);
			};

			auto chi2_landscape_pane = gui::vnotebook(
				view,
				gui::deck(
					make_image_pane(settings::general::output + plot_names[0].first + ".png"),
					make_image_pane(settings::general::output + plot_names[1].first + ".png"),
					make_image_pane(settings::general::output + plot_names[2].first + ".png"),
					make_image_pane(settings::general::output + plot_names[3].first + ".png"),
					make_image_pane(settings::general::output + plot_names[4].first + ".png")
				),
				gui::tab(plot_names[0].second),
				gui::tab(plot_names[1].second),
				gui::tab(plot_names[2].second),
				gui::tab(plot_names[3].second),
				gui::tab(plot_names[4].second)
			);

			auto chi2_pane = gui::vnotebook(
				view,
				gui::deck(
					make_image_pane(settings::general::output + plot_names[5].first + ".png"),
					make_image_pane(settings::general::output + plot_names[6].first + ".png")
				),
				gui::tab("log"),
				gui::tab("log-log")
			);

			auto image_viewer_layout = gui::margin(
				{10, 10, 10, 10},
				gui::align_center_middle(
					gui::vnotebook(
						view,
						gui::deck(
							chi2_pane,
							chi2_landscape_pane
						),
						gui::tab("Scattering profile"),
						gui::tab("χ² landscape")
					)
				)
			);

			deck.push_back(gui::share(image_viewer_layout));
			deck.select(2);
			view.refresh();
		});
	};

	return link(deck);
}

int main(int, char*[]) {
    std::ios_base::sync_with_stdio(false);
	settings::axes::qmin = 0;
	settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;
    settings::em::mass_axis = true;
    settings::em::hydrate = true;
    settings::fit::verbose = true;
    settings::em::alpha_levels = {1, 10};
    settings::hist::weighted_bins = true;

	resources::generate_resource_file();
	auto logo_path = resources::generate_logo_file();

	gui::app app("EM fitter");
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
				gui::label(std::string(constants::version))
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
				gui::scale_element(0.15, gui::image(abs_path(logo_path).c_str()))
			)
		)
	);
	auto settings = gui::vtile(
		gui::margin_top(
			10,
			gui::label("Input & output")
		),
		io_menu(view),
		gui::hgrid(
			{0.5, 1.0},
			gui::vtile(
				gui::margin_top(
					10,
					gui::label("q-range")
				),
				q_slider(view)
			),
			gui::vtile(
				gui::margin_top(
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
