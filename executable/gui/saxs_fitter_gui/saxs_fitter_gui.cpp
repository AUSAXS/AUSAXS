// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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

#include <io/File.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <fitter/FitReporter.h>
#include <fitter/SmartFitter.h>
#include <utility/Logging.h>
#include <settings/All.h>
#include <plots/All.h>

#include <gui/helper.h>
#include <gui/logo.h>
#include <gui/resources.h>

#include <bitset>
#include <thread>
#include <string_view>

using namespace ausaxs;
namespace gui = cycfi::elements;

auto plot_names = std::vector<std::pair<std::string, std::string>>{
	{"log", "single-logarithmic plot"},
	{"loglog", "double-logarithmic plot"},
	{"p(r)", "distance histogram"},
	{"profiles", "partial profiles"}
};

auto make_start_button(gui::view& view) {
	static auto start_button = gui::button("Start");
	start_button->set_body_color(ColorManager::get_color_success());

	auto start_button_layout = gui::margin(
		{10, 100, 10, 100},
		gui::align_center_middle(
			gui::hsize(
				200,
				link(start_button)
			)
		)
	);

	static auto deck = gui::deck(
		start_button_layout,
		start_button_layout
	);

	static std::thread worker;
	start_button.on_click = [&view] (bool) {
		// ensure the worker is ready to be assigned a job
		if (worker.joinable()) {
			worker.join();
		}

		if (!setup::saxs_dataset || !io::File(::settings::pdb_file).exists()) {
			ausaxs::console::print_warning("no saxs data or pdb file was provided");
			start_button->set_body_color(ColorManager::get_color_fail());
			start_button->set_text("missing input");
			view.refresh(start_button);
			worker = std::thread([&view] () {
				std::this_thread::sleep_for(std::chrono::milliseconds(2000));
				start_button->set_body_color(ColorManager::get_color_success());
				start_button->set_text("start");
				view.refresh(start_button);
			});
			return;
		}

		start_button->set_body_color(ColorManager::get_color_accent());
		start_button->set_text("Working...");

		// use a worker thread to avoid locking the gui
		worker = std::thread([&view] () {
			ausaxs::settings::axes::qmin = qslider_axis_transform(qslider->value_first());
			ausaxs::settings::axes::qmax = qslider_axis_transform(qslider->value_second());
			logging::log(std::string("Starting SAXS fitting with: ") +
				"\n\tqmin = " + std::to_string(ausaxs::settings::axes::qmin) +
				"\n\tqmax = " + std::to_string(ausaxs::settings::axes::qmax) +
				"\n\tpdb file = " + ::settings::pdb_file +
				"\n\tsaxs file = " + ::settings::saxs_file
			);
			setup::pdb = std::make_unique<data::Molecule>(::settings::pdb_file);
			setup::pdb->generate_new_hydration();

			fitter::SmartFitter fitter({::settings::saxs_file}, setup::pdb->get_histogram());
			auto result = fitter.fit();

			fitter::FitReporter::report(result.get());
			fitter::FitReporter::save(result.get(), ausaxs::settings::general::output + "report.txt");

			plots::PlotDistance::quick_plot(fitter.get_model(), ausaxs::settings::general::output + "p(r)." + ausaxs::settings::plots::format);
			plots::PlotProfiles::quick_plot(fitter.get_model(), ausaxs::settings::general::output + "profiles." + ausaxs::settings::plots::format);
			result->curves.save(ausaxs::settings::general::output + "ausaxs.fit", "chi2=" + std::to_string(result->fval/result->dof) + " dof=" + std::to_string(result->dof));

			setup::pdb->save(ausaxs::settings::general::output + "model.pdb");
			perform_plot(ausaxs::settings::general::output);

			auto make_image_pane = [] (const io::File& path) {
				if (!path.exists()) {
					throw except::io_error("File " + path.absolute_path() + " does not exist");
				}
				logging::log("Loading image " + path.absolute_path());
				return gui::image(path.absolute_path().c_str(), 0.13);
			};

			auto main_pane = gui::vnotebook(
				view,
				gui::deck(
					make_image_pane(ausaxs::settings::general::output + plot_names[0].first + ".png"),
					make_image_pane(ausaxs::settings::general::output + plot_names[1].first + ".png"),
					make_image_pane(ausaxs::settings::general::output + plot_names[2].first + ".png"),
					make_image_pane(ausaxs::settings::general::output + plot_names[3].first + ".png")
				),
				gui::tab(plot_names[0].second),
				gui::tab(plot_names[1].second),
				gui::tab(plot_names[2].second),
				gui::tab(plot_names[3].second)
			);

			auto image_viewer_layout = 	gui::margin(
				{10, 10, 10, 10},
				gui::vtile(
					main_pane,
					gui::margin_top(
						10,
						gui::align_center_middle(
							gui::hsize(
								200,
								link(start_button)
							)
						)
					)
				)
			);
			logging::log("Created image viewer layout");

			start_button->set_body_color(ColorManager::get_color_success());
			start_button->set_text("Start");

			logging::log("Switching to image viewer layout");
			deck[1] = gui::share(image_viewer_layout);
			deck.select(1);
			view.refresh();
			logging::log("Finished processing SAXS data and updating GUI");
		});
	};

	return link(deck);
}

auto io_menu(gui::view& view) {
	static auto saxs_box_bg = ColorManager::new_background_color();
	static auto pdb_box_bg = ColorManager::new_background_color();
	static auto output_box_bg = ColorManager::new_background_color();

	static auto saxs_box = gui::input_box("SAXS path");
	static auto pdb_box = gui::input_box("Structure path");
	static auto output_box = gui::input_box("Output path");
	ColorManager::manage_input_box(saxs_box.second);
	ColorManager::manage_input_box(pdb_box.second);
	ColorManager::manage_input_box(output_box.second);

	static bool default_output = true;
	output_box.second->set_text("output/saxs_fitter");
	static bool pdb_ok = false, saxs_ok = false;

	pdb_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			pdb_box_bg.get()._color = ColorManager::get_color_accent();
		} else if (text.empty()) {
			pdb_box_bg.get()._color = ColorManager::get_color_background();
		}

		if (pdb_ok) {
			pdb_ok = false;
			pdb_box_bg.get() = ColorManager::get_color_accent();
		}

		static unsigned int last_size = 0;
		auto fill = autocomplete(text, last_size, [] (const io::File& p) {return constants::filetypes::structure.check(p);});
		if (!fill.first.empty()) {pdb_box.second->set_text(fill.first);}
		if (fill.second) {pdb_box.second->on_enter(fill.first);}
	};

	pdb_box.second->on_enter = [] (std::string_view text) -> bool {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::structure.check(file)) {
			console::print_warning("invalid pdb file " + file.path());
			pdb_box_bg.get() = ColorManager::get_color_fail();
			pdb_ok = false;
			return true;
		}

		// check if we can use a relative path instead of absolute
		if (auto curpath = std::filesystem::current_path().string(); file.path().find(curpath) != std::string::npos) {
			file = std::filesystem::relative(file.path(), curpath).string();
		}
		pdb_box.second->set_text(file.path());

		::settings::pdb_file = file.path();
		logging::log("PDB file was set to " + ::settings::pdb_file);
		pdb_box_bg.get() = ColorManager::get_color_success();
		pdb_ok = true;

		if (!saxs_ok) {
			if (20 < std::distance(std::filesystem::directory_iterator(file.directory().path()), std::filesystem::directory_iterator{})) {return true;}
			for (auto& p : std::filesystem::directory_iterator(file.directory().path())) {
				io::File tmp(p.path().string());
				if (constants::filetypes::saxs_data.check(tmp)) {
					::settings::saxs_file = tmp.path();
					saxs_box.second->set_text(tmp.path());
					saxs_box.second->on_enter(tmp.path());
				}
			}
		}

		if (saxs_ok && default_output) {
			std::string path = "output/saxs_fitter/" + io::File(::settings::pdb_file).stem() + "/" + io::File(::settings::saxs_file).stem();
			output_box.second->set_text(path);
			output_box.second->on_enter(path);
		}
		return true;
	};

	saxs_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			saxs_box_bg.get() = ColorManager::get_color_accent();
		} else if (text.empty()) {
			saxs_box_bg.get() = ColorManager::get_color_background();
		}

		if (saxs_ok) {
			saxs_ok = false;
			saxs_box_bg.get() = ColorManager::get_color_accent();
		}

		static unsigned int last_size = 0;
		auto fill = autocomplete(text, last_size, [] (const io::File& p) {return constants::filetypes::saxs_data.check(p);});
		if (!fill.first.empty()) {saxs_box.second->set_text(fill.first);}
		if (fill.second) {saxs_box.second->on_enter(fill.first);}
	};

	saxs_box.second->on_enter = [] (std::string_view text) -> bool {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::saxs_data.check(file)) {
			console::print_warning("invalid saxs file " + file.path());
			saxs_box_bg.get() = ColorManager::get_color_fail();
			saxs_ok = false;
			setup::saxs_dataset = nullptr;
			return true;
		}

		// check if we can use a relative path instead of absolute
		if (auto curpath = std::filesystem::current_path().string(); file.path().find(curpath) != std::string::npos) {
			file = std::filesystem::relative(file.path(), curpath).string();
		}
		saxs_box.second->set_text(file.path());

		logging::log("SAXS file was set to " + ::settings::saxs_file);
		::settings::saxs_file = file.path();
		saxs_box_bg.get() = ColorManager::get_color_success();
		setup::saxs_dataset = std::make_unique<SimpleDataset>(::settings::saxs_file);
		saxs_ok = true;

		if (pdb_ok) {
		 	if (default_output || output_box.second->get_text().empty()) {
				std::string path = "output/saxs_fitter/" + io::File(::settings::pdb_file).stem() + "/" + io::File(::settings::saxs_file).stem();
				output_box.second->set_text(path);
				output_box.second->on_enter(path);
			}
		}
		qslider_update_range_from_file();
		return true;
	};

	output_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			output_box_bg.get() = ColorManager::get_color_accent();
		} else if (text.empty()) {
			output_box_bg.get() = ColorManager::get_color_background();
		}
		default_output = false;
	};

	output_box.second->on_enter = [&view] (std::string_view text) -> bool {
		ausaxs::settings::general::output = text;
		if (ausaxs::settings::general::output.back() != '/') {
			ausaxs::settings::general::output += "/";
			view.refresh(output_box.first);
		}
		console::print_text("output path was set to " + ausaxs::settings::general::output);
		return true;
	};

	auto map_box_field = make_file_dialog_button(pdb_box, pdb_box_bg, {"PDB file", "pdb"});
	auto saxs_box_field = make_file_dialog_button(saxs_box, saxs_box_bg, {"SAXS data", "dat,scat"});
	auto output_box_field = make_folder_dialog_button(output_box, output_box_bg);

	return gui::htile(
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
	);
}

auto selection_menu_settings(gui::view&) {
	// we use a deck composite to avoid circular dependencies
	static auto deck = gui::deck_composite();

	std::vector<std::pair<std::string, ausaxs::settings::hydrate::HydrationStrategy>> options1 {
		{"1. Radial", ausaxs::settings::hydrate::HydrationStrategy::RadialStrategy},
		{"2. None", ausaxs::settings::hydrate::HydrationStrategy::NoStrategy}
	};
	static auto hydration_model = gui::selection_menu(
		[options1] (std::string_view selection) {
			for (auto& option : options1) {
				if (option.first == selection) {
					ausaxs::settings::hydrate::hydration_strategy = option.second;
				}
			}
		}, 
		{
			options1[0].first,
			options1[1].first
		}
	);
	ColorManager::manage_text([] (gui::color color) {hydration_model.second->font_color(color);});

	std::vector<std::pair<std::string, ausaxs::settings::exv::ExvMethod>> options2 {
		{"1. Simple", ausaxs::settings::exv::ExvMethod::Simple},
		{"2. Fraser", ausaxs::settings::exv::ExvMethod::Fraser},
		{"3. Grid", ausaxs::settings::exv::ExvMethod::GridSurface}
	};

	static auto excluded_volume_model = gui::selection_menu(
		[options2] (std::string_view selection) {
			for (auto& option : options2) {
				if (option.first == selection) {
					ausaxs::settings::exv::exv_method = option.second;
				}
			}
			switch (ausaxs::settings::exv::exv_method) {
				case ausaxs::settings::exv::ExvMethod::Simple:
					ausaxs::settings::fit::fit_excluded_volume = false;
					deck.select(0);
					break;
				default:
					deck.select(1);
					break;
			}
		}, 
		{
			options2[0].first,
			options2[1].first,
			options2[2].first
		}
	);
	ColorManager::manage_text([] (gui::color color) {excluded_volume_model.second->font_color(color);});

	static auto fit_excluded_volume_button = gui::check_box("fit excluded volume");
	fit_excluded_volume_button.on_click = [] (bool value) {
		ausaxs::settings::fit::fit_excluded_volume = value;
	};

	static auto fit_solvent_density_button = gui::check_box("fit solvent density");
	fit_solvent_density_button.on_click = [] (bool value) {
		ausaxs::settings::fit::fit_solvent_density = value;
	};

	static auto hydration_text = gui::label("Hydration model")
		.font_color(ColorManager::get_text_color())
		.font_size(18);
	ColorManager::manage_text([] (gui::color color) {hydration_text.set_font_color(color);});
	static auto exv_text = gui::label("Excluded volume model")
		.font_color(ColorManager::get_text_color())
		.font_size(18);
	ColorManager::manage_text([] (gui::color color) {exv_text.set_font_color(color);});
	
	auto exv_fit_support_layout = gui::margin_left_right(
		{10, 10},
		gui::htile(
			gui::vtile(
				gui::margin_bottom(
					10,
					link(hydration_text)
				),
				link(hydration_model.first)
			),
			gui::hspace(50),
			gui::vtile(
				gui::margin_bottom(
					10,
					link(exv_text)
				),
				link(excluded_volume_model.first)
			),
			gui::hspace(50),
			gui::vtile(
				gui::margin(
					{10, 10, 0, 10},
					gui::fixed_size(
						{200, 100},
						link(fit_excluded_volume_button)
					)
				),
				gui::margin(
					{10, 10, 0, 10},
					gui::fixed_size(
						{200, 100},
						link(fit_solvent_density_button)
					)
				)
			)
		)
	);

	auto no_exv_fit_support_layout = gui::margin_left_right(
		{10, 10},
		gui::htile(
			gui::vtile(
				gui::margin_bottom(
					10,
					link(hydration_text)
				),
				link(hydration_model.first)
			),
			gui::hspace(50),
			gui::vtile(
				gui::margin_bottom(
					10,
					link(exv_text)
				),
				link(excluded_volume_model.first)
			)
		)
	);

	deck.push_back(share(no_exv_fit_support_layout));
	deck.push_back(share(exv_fit_support_layout));
	deck.select(0);

	return link(deck);
}

// toggle light/dark mode
auto toggle_mode_button(gui::view& view) {
	static auto button = gui::button("light mode");
	button.on_click = [&view] (bool) {
		button->set_text(ColorManager::dark_mode ? "light mode" : "dark mode");
		ColorManager::switch_mode();
		view.refresh();
	};
	return link(button);
}

int main(int, char*[]) {
    std::ios_base::sync_with_stdio(false);

	logging::start("saxs_fitter_gui");
	console::print_info("Starting saxs_fitter_gui");

	gui::app app("AUSAXS saxs fitter");
	gui::window win(app.name(), std::bitset<4>{"1111"}.to_ulong(), {50, 50, 1024+50, 768+50});
	win.on_close = [&app]() {app.stop(); exit(0);};

	resources::generate_resource_file();
	auto logo_path = resources::generate_logo_file();

	auto title = gui::label("SAXS fitter")
		.font_size(50)
		.font_color(ColorManager::get_text_color());
	ColorManager::manage_text([&title] (gui::color color) {title.set_font_color(color);});

	gui::view view(win);
	auto header = gui::layer(
		gui::align_center_top(
			link(title)
		),
		gui::align_right_top(
			gui::margin(
				{0, 10, 10, 10},
				gui::fixed_size(
					{ 50, 50 },
					toggle_mode_button(view)
				)
			)
		)
	);
	auto content = gui::margin(
		{10, 10, 10, 10},
		gui::vtile(
			io_menu(view),
			gui::htile(
				q_slider(view),
				selection_menu_settings(view)
			),
			gui::align_center_middle(
				make_start_button(view)
			)
		)
	);
	auto authors = gui::label("Kristian Lytje & Jan Skov Pedersen")
		.font_color(ColorManager::get_text_color());
	ColorManager::manage_text([&authors] (gui::color color) {authors.set_font_color(color);});
	auto version = gui::label(std::string(constants::version))
		.font_color(ColorManager::get_text_color());
	ColorManager::manage_text([&version] (gui::color color) {version.set_font_color(color);});
	auto footer = gui::margin(
		{10, 10, 10, 10},
		gui::hgrid(
			{0.33, 0.66, 1},
			gui::align_left_bottom(
				link(version)
			),
			gui::align_center_bottom(
				link(authors)
			),
			gui::align_right_bottom(
				gui::scale_element(0.15, gui::image(logo_path.absolute_path().c_str()))
			)
		)
	);

	view.content(
		gui::vtile(
			header,
			content,
			footer
		),
		link(background)
	);

	app.run();
	return 0;
}
