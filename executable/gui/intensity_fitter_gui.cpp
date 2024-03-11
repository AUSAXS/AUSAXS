#include "fitter/ExcludedVolumeFitter.h"
#include <elements.hpp>
#include <nfd.hpp>
#include <include.h>

#include <io/File.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <fitter/HydrationFitter.h>
#include <fitter/FitReporter.h>
#include <fitter/Fit.h>
#include <settings/All.h>
#include <plots/All.h>

#include <bitset>

namespace gui = cycfi::elements;

auto plot_names = std::vector<std::pair<std::string, std::string>>{
	{"log", "single-logarithmic plot"},
	{"loglog", "double-logarithmic plot"},
	{"p(r)", "distance histogram"},
	{"profiles", "partial profiles"}
};

auto make_start_button(gui::view& view) {
	static auto start_button = gui::button("start");

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

	static std::thread worker;
	start_button.on_click = [&view] (bool) {
		if (!setup::saxs_dataset || !io::File(settings::pdb_file).exists()) {
			std::cout << "no saxs data or pdb file was provided" << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			return;
		}
		setup::pdb = std::make_unique<data::Molecule>(settings::pdb_file);
		view.refresh();
		worker = std::thread([&view] () {
			bool fit_excluded_volume = false; //!

			std::shared_ptr<fitter::HydrationFitter> fitter;
			if (fit_excluded_volume) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(settings::saxs_file, setup::pdb->get_histogram());}
			else {fitter = std::make_shared<fitter::HydrationFitter>(settings::saxs_file, setup::pdb->get_histogram());}
			std::shared_ptr<fitter::Fit> result = fitter->fit();

			fitter::FitReporter::report(result.get());
			fitter::FitReporter::save(result.get(), settings::general::output + "report.txt");

			plots::PlotDistance::quick_plot(fitter->get_scattering_hist(), settings::general::output + "p(r)." + settings::plots::format);
			plots::PlotProfiles::quick_plot(fitter->get_scattering_hist(), settings::general::output + "profiles." + settings::plots::format);

			fitter->get_model_dataset().save(settings::general::output + "fit.fit");
			fitter->get_dataset().save(settings::general::output + io::File(settings::saxs_file).stem() + ".scat");

			setup::pdb->save(settings::general::output + "model.pdb");
			perform_plot(settings::general::output);

			auto make_image_pane = [] (const io::File& path) {
				return gui::image(std::filesystem::current_path().string() + "/" + path.path().c_str(), 0.15);
			};

			auto main_pane = gui::vnotebook(
				view,
				gui::deck(
					make_image_pane(settings::general::output + plot_names[0].first + ".png"),
					make_image_pane(settings::general::output + plot_names[1].first + ".png"),
					make_image_pane(settings::general::output + plot_names[2].first + ".png"),
					make_image_pane(settings::general::output + plot_names[3].first + ".png")
				),
				gui::tab(plot_names[0].second),
				gui::tab(plot_names[1].second),
				gui::tab(plot_names[2].second),
				gui::tab(plot_names[3].second)
			);

			auto image_viewer_layout = gui::margin(
				{10, 10, 10, 10},
				main_pane
			);

			deck.push_back(gui::share(image_viewer_layout));
			deck.select(1);
		});
	};

	return link(deck);
}

auto io_menu(gui::view& view) {
	static auto saxs_box_bg = ColorManager::new_background_color();
	static auto pdb_box_bg = ColorManager::new_background_color();
	static auto output_box_bg = ColorManager::new_background_color();

	static auto saxs_box = gui::input_box("saxs path");
	static auto pdb_box = gui::input_box("pdb path");
	static auto output_box = gui::input_box("output path");
	saxs_box.second->set_color(ColorManager::get_text_color());
	pdb_box.second->set_color(ColorManager::get_text_color());
	output_box.second->set_color(ColorManager::get_text_color());

	static bool default_output = true;
	output_box.second->set_text("output/saxs_fitter");
	auto ref = output_box.second;
	static bool pdb_ok = false, saxs_ok = false;

	pdb_box.second->on_text = [&view] (std::string_view text) {
		static unsigned int last_size = 0;
		if (text.size() == 1) {
			pdb_box_bg.get()._color = ColorManager::get_color_accent();
		} else if (text.empty()) {
			pdb_box_bg.get()._color = ColorManager::get_color_background();
		}

		if (pdb_ok) {
			pdb_ok = false;
			pdb_box_bg.get() = ColorManager::get_color_accent();
		}

		// prevent autocompletion when deleting text
		if (text.size() < last_size) {
			last_size = text.size();
			return;
		}
		last_size = text.size();

		// only autocomplete if the last character is a '/' and there are less than 20 matches
		if (text.back() != '/') {return;}
		if (!std::filesystem::is_directory(text)) {return;}
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
			// settings::map_file = matches.front();
			pdb_box.second->set_text(matches.front());
			pdb_box.second->on_enter(matches.front());
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
			pdb_box.second->set_text(prefix);
		}
	};

	pdb_box.second->on_enter = [&view] (std::string_view text) {
		io::File file = io::File(std::string(text));
		if (!constants::filetypes::structure.validate(file)) {
			std::cout << "invalid pdb file " << file.path() << std::endl;
			pdb_box_bg.get() = ColorManager::get_color_fail();
			pdb_ok = false;
			return;
		}

		settings::pdb_file = file.path();
		std::cout << "pdb file was set to " << settings::pdb_file << std::endl;
		pdb_box_bg.get() = ColorManager::get_color_success();
		pdb_ok = true;

		if (!saxs_ok) {
			if (20 < std::distance(std::filesystem::directory_iterator(file.directory().path()), std::filesystem::directory_iterator{})) {return;}
			for (auto& p : std::filesystem::directory_iterator(file.directory().path())) {
				io::File tmp(p.path().string());
				if (constants::filetypes::saxs_data.validate(tmp)) {
					settings::saxs_file = tmp.path();
					saxs_box.second->set_text(tmp.path());
					saxs_box.second->on_enter(tmp.path());
				}
			}
		}

		if (saxs_ok && default_output) {
			std::string path = "output/saxs_fitter/" + io::File(settings::pdb_file).stem() + "/" + io::File(settings::saxs_file).stem();
			output_box.second->set_text(path);
			output_box.second->on_enter(path);
		}
	};

	saxs_box.second->on_text = [&view] (std::string_view text) {
		static unsigned int last_size = 0;
		if (text.size() == 1) {
			saxs_box_bg.get() = ColorManager::get_color_accent();
		} else if (text.empty()) {
			saxs_box_bg.get() = ColorManager::get_color_background();
		}

		if (saxs_ok) {
			saxs_ok = false;
			saxs_box_bg.get() = ColorManager::get_color_accent();
		}

		// prevent autocompletion when deleting text
		if (text.size() < last_size) {
			last_size = text.size();
			return;
		}
		last_size = text.size();

		// only autocomplete if the last character is a '/' and there are less than 20 matches
		if (text.back() != '/') {return;}
		if (!std::filesystem::is_directory(text)) {return;}
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
			saxs_box_bg.get() = ColorManager::get_color_fail();
			saxs_ok = false;
			setup::saxs_dataset = nullptr;
			return;
		}

		settings::saxs_file = file.path();
		std::cout << "saxs file was set to " << settings::saxs_file << std::endl;
		saxs_box_bg.get() = ColorManager::get_color_success();
		saxs_ok = true;
		setup::saxs_dataset = std::make_unique<SimpleDataset>(settings::saxs_file);

		if (pdb_ok) {
		 	if (default_output) {
				std::string path = "output/saxs_fitter/" + io::File(settings::pdb_file).stem() + "/" + io::File(settings::saxs_file).stem();
				output_box.second->set_text(path);
				output_box.second->on_enter(path);
			}
		}
	};

	output_box.second->on_text = [] (std::string_view text) {
		if (text.size() == 1) {
			output_box_bg.get() = ColorManager::get_color_accent();
		} else if (text.empty()) {
			output_box_bg.get() = ColorManager::get_color_background();
		}
		default_output = false;
	};

	output_box.second->on_enter = [&view] (std::string_view text) {
		settings::general::output = text;
		if (settings::general::output.back() != '/') {
			settings::general::output += "/";
			view.refresh(output_box.first);
		}
		std::cout << "output path was set to " << settings::general::output << std::endl;
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

// toggle light/dark mode
auto toggle_mode_button(gui::view& view) {
	static auto button = gui::button("light mode");
	button.on_click = [&view] (bool checked) {
		ColorManager::switch_mode();
		button->set_text(ColorManager::dark_mode ? "light mode" : "dark mode");
		view.refresh();
	};
	return link(button);
}

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
	gui::app app(argc, argv, "AUSAXS intensity fitter", "com.cycfi.ausaxs-intensity-fitter");
	gui::window win(app.name(), std::bitset<4>{"1111"}.to_ulong());
	win.on_close = [&app]() {app.stop();};

	io::File logo_path = "logo.png";

	gui::view view(win);
	auto header = gui::layer(
		gui::align_center_top(
			gui::label("Intensity fitter")
				.font_size(50)
				.font_color(ColorManager::get_text_color())
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
			q_slider(view),
			gui::align_center_middle(
				make_start_button(view)
			)
		)
	);
	auto footer = gui::margin(
		{10, 10, 10, 10},
		gui::hgrid(
			{0.33, 0.66, 1},
			gui::align_left_bottom(
				gui::label(std::string(constants::version)).font_color(ColorManager::get_text_color())
			),
			gui::align_center_bottom(
				gui::label("Kristian Lytje & Jan Skov Pedersen").font_color(ColorManager::get_text_color())
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
