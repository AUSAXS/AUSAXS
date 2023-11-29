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

using namespace cycfi::elements;

// Main window background color
auto constexpr bkd_color_accent = rgba(55, 55, 57, 255);
auto constexpr bkd_color = rgba(35, 35, 37, 255);
constexpr auto bgreen   = colors::green.level(0.7).opacity(0.4);
auto background = box(bkd_color);

// The different image types which can be plotted
std::shared_ptr<cycfi::elements::image> image_intensity, image_distance, image_intensity_fit, image_residuals; 

// The currently shown image
std::shared_ptr<cycfi::elements::image> current_image;
int current_image_selection = 1;

// Sensible defaults
std::string file_in = "data/2epe.pdb", file_out = "figures/";

auto io_menu() {
   static float const grid[] = { 0.32, 1.0 };

   auto my_label = [=](auto text)
   {
      return right_margin(10, label(text).text_align(canvas::right));
   };

   auto my_input = [=] (auto caption, auto input)
   {
      return bottom_margin(10, hgrid(grid, my_label(caption), input));
   };

   // This is an example on how to add an on_text callback:
   auto in = input_box("Input path");
   in.second->set_text("2epe");
   in.second->on_text = [input = in.second.get()] (std::string_view text)
      {
         file_in = "data/" + std::string(text) + ".pdb";
         std::cout << "Input file is currently " << file_in << std::endl;
      };

   auto out = input_box("Output path");
   out.second->set_text("figures/");
   out.second->on_text = [input = out.second.get()] (std::string_view text)
      {
         file_out = text;
         std::cout << "Output file is currently " << file_out << std::endl;
      };

   return htile(
               my_input("Input path", in.first), 
               left_margin(20, my_input("Output path", out.first))
          );
}

auto placement_strategy_menu()
{
   return selection_menu(
      [] (std::string_view select)
      {
         std::cout << "Selected " << select << std::endl;
         if (select == "Radial strategy") {
            settings::grid::placement_strategy = settings::grid::PlacementStrategy::RadialStrategy;
         } else if (select == "Axes strategy") {
            settings::grid::placement_strategy = settings::grid::PlacementStrategy::AxesStrategy;
         } else if (select == "Jan strategy") {
            settings::grid::placement_strategy = settings::grid::PlacementStrategy::JanStrategy;
         }
      },
      {
         "Radial strategy",
         "Axes strategy",
         "Jan strategy"
      }

   ).first; // We return only the first, the menu. the second is a shared pointer to the label.
}

auto widths_menu() 
{
   static float const grid[] = {0.4, 1.0};

   auto my_label = [=](auto text)
   {
      return right_margin(10, label(text).text_align(canvas::right));
   };

   auto my_input = [=] (auto caption, auto input)
   {
      return bottom_margin(10, hgrid(grid, my_label(caption), input));
   };

   auto bin_width = input_box("Bin width");
   bin_width.second->set_text(std::to_string(constants::axes::q_axis.bins));
   bin_width.second->on_text = [input = bin_width.second.get()] (std::string_view text)
      {
         try {
         } catch (const std::exception&) {}
      };

   auto grid_width = input_box("Grid width");
   grid_width.second->set_text(std::to_string(settings::grid::width));
   grid_width.second->on_text = [input = grid_width.second.get()] (std::string_view text)
      {
         try {
            settings::grid::width = std::stod(std::string(text));
         } catch (const std::exception&) {}
      };

   return htile(
               my_input("Bin width", bin_width.first), 
               left_margin(20, my_input("Grid width", grid_width.first))
          );
}

auto toggle_menu(view& _view) {
   auto check_center = check_box("Center molecule");
   auto check_effective_charge = check_box("Use effective charge");
   check_center.value(true);
   check_effective_charge.value(true);

   check_center.on_click = [&_view] (bool pressed) mutable
   {
      settings::molecule::center = pressed;
   };

   check_effective_charge.on_click = [&_view] (bool pressed) mutable
   {
      settings::molecule::use_effective_charge = pressed;
      if (pressed) {
         std::cout << "Enabled effective charge" << std::endl;
      } else {
         std::cout << "Disabled effective charge" << std::endl;
      }
   };

   return group("Check boxes",
            margin({ 10, 10, 20, 20 },
               top_margin(25,
                  htile(
                     top_margin(10, align_left(check_center)),
                     top_margin(10, align_left(check_effective_charge))
                  )
               )
            )
   );
}

auto radii_menu() 
{
   static float const grid[] = {0.4, 1.0};

   auto my_label = [=](auto text)
   {
      return right_margin(10, label(text).text_align(canvas::right));
   };

   auto my_input = [=] (auto caption, auto input)
   {
      return bottom_margin(10, hgrid(grid, my_label(caption), input));
   };

   auto ra = input_box("Excluded volume radius");
   ra.second->set_text(std::to_string(settings::grid::rvol));
   ra.second->on_text = [input = ra.second.get()] (std::string_view text)
      {
         try {
            settings::grid::rvol = std::stod(std::string(text));
         } catch (const std::exception&) {}
      };

   auto rh = input_box("Water radius");
   rh.second->set_text(std::to_string(settings::grid::rvol));
   rh.second->on_text = [input = rh.second.get()] (std::string_view text)
      {
         try {
            settings::grid::rvol = std::stod(std::string(text));
         } catch (const std::exception&) {}
      };

   return htile(
               my_input("Atomic radius", ra.first), 
               left_margin(20, my_input("Water radius", rh.first))
          );
}

auto image_control(view& _view) {
   auto radio_distance = radio_button("Distance histogram");
   auto radio_intensity = radio_button("Scattering intensity");
   auto radio_intensity_fit = radio_button("Scattering intensity fit");
   auto radio_residuals = radio_button("Fit residuals");
   radio_distance.select(true);

   radio_distance.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          *current_image = *image_distance;
          current_image_selection = 1;
          _view.refresh();
      }
   };

   radio_intensity.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          *current_image = *image_intensity;
          current_image_selection = 2;
          _view.refresh();
      }
   };

   radio_intensity_fit.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          *current_image = *image_intensity_fit;
          current_image_selection = 3;
          _view.refresh();
      }
   };

   radio_residuals.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          *current_image = *image_residuals;
          current_image_selection = 4;
          _view.refresh();
      }
   };

   return group("Change the shown figure",
            margin({10, 10, 20, 20},
               top_margin(25,
                  htile
                  (
                     top_margin(10, align_left(radio_distance)),
                     top_margin(10, align_left(radio_intensity)),
                     top_margin(10, align_left(radio_intensity_fit)),
                     top_margin(10, align_left(radio_residuals))
                  )
               )
            )
         );
}

auto generate_button(view& _view) {
   auto lbutton = button("Generate figures", 1.0, bgreen);
   lbutton.on_click = [&_view] (bool pressed) mutable
      {
         if (pressed) {
            data::Molecule protein(file_in);
            protein.generate_new_hydration();
            auto d = protein.get_histogram();
            auto p = d->debye_transform();

            // Distance plot
            std::string path = "temp/gui/distance.png";
            plots::PlotDistance d_plot(d.get(), path);
            d_plot.save(path); 
            image_distance = share(image("temp/gui/distance.png"));

            // Debye scattering intensity plot
            plots::PlotIntensity i_plot(p, plots::PlotOptions());
            path = "temp/gui/log.png";
            i_plot.save(path);
            image_intensity = share(image("temp/gui/log.png"));

            // Debye scattering intensity fit
            std::string measurement_data = file_in.substr(0, file_in.find_last_of('.')) + ".rsr";
            fitter::HydrationFitter fitter(measurement_data, std::move(d));
            auto result = fitter.fit();

            plots::PlotIntensityFit plot_f(fitter);
            path = "temp/gui/log.png";
            plot_f.save(path);
            image_intensity_fit = share(image("temp/gui/log.png"));

            plots::PlotIntensityFitResiduals plot_r(fitter);
            path = "temp/gui/log.png";
            plot_r.save(path);
            image_residuals = share(image("temp/gui/log.png"));

            if (current_image_selection == 1) {*current_image = *image_distance;} 
            else if (current_image_selection == 2) {*current_image = *image_intensity;} 
            else if (current_image_selection == 3) {*current_image = *image_intensity_fit;} 
            else if (current_image_selection == 4) {*current_image = *image_residuals;} 

            _view.refresh();
         }
      };
   return lbutton;
}

auto make_controls(view& _view) {
   image_intensity     = share(image(std::string(std::filesystem::current_path().string() + "/temp/gui/log.png").c_str()));
   image_intensity_fit = share(image(std::string(std::filesystem::current_path().string() + "/temp/gui/log.png").c_str()));
   image_residuals     = share(image(std::string(std::filesystem::current_path().string() + "/temp/gui/log.png").c_str()));
   image_distance      = share(image(std::string(std::filesystem::current_path().string() + "/temp/gui/distance.png").c_str()));
   current_image       = share(image(std::string(std::filesystem::current_path().string() + "/temp/gui/distance.png").c_str())); //! Must be assigned this way - otherwise compiler optimizations may remove later assignments
   auto image_box = layer
                        (
                           margin
                           (
                              {10, 10, 10, 10},
                              link(*current_image)
                           ),
                           box(bkd_color_accent)
                        );

    return vtile
            (
               margin({10, 5, 5, 10}, io_menu()),
               margin({10, 5, 5, 10}, placement_strategy_menu()),
               margin({10, 5, 5, 10}, widths_menu()),
               margin({10, 5, 5, 10}, radii_menu()),
               margin({10, 5, 5, 10}, toggle_menu(_view)),
               htile
               (
                  margin({10, 20, 5, 10}, image_control(_view))
               ),
               left_margin(10, generate_button(_view)),
               image_box
            );
}

int main(int argc, char* argv[]) {
   app _app(argc, argv, "SAXS", "com.saxs.gui");
   window _win(_app.name(), 15, rect{20, 20, 1000, 1500});
   _win.on_close = [&_app]() { _app.stop(); };

   view _view(_win);

   _view.content
   (
      make_controls(_view),
      background
   );

   _app.run();
   return 0;
}
