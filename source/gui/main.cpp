/*=============================================================================
   Copyright (c) 2016-2020 Joel de Guzman

   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/
#include <elements.hpp>
#include <algorithm>
#include <random>
#include <iostream>

#include "data/Protein.h"
#include "plots/PlotDistance.h"
#include "plots/PlotIntensity.h"
#include "settings.h"

using namespace cycfi::elements;
using std::cout, std::endl; 

// Main window background color
auto constexpr bkd_color_accent = rgba(55, 55, 57, 255);
auto constexpr bkd_color = rgba(35, 35, 37, 255);
constexpr auto bgreen   = colors::green.level(0.7).opacity(0.4);
auto background = box(bkd_color);

// The different image types which can be plotted
std::shared_ptr<cycfi::elements::image> image_intensity, image_distance; 

// The currently shown image
std::shared_ptr<cycfi::elements::image> current_image;
int current_image_selection = 1;

// Sensible defaults
string file_in = "data/2epe.pdb", file_out = "figures/";

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
   in.second->on_text =
      [input = in.second.get()](std::string_view text)
      {
         file_in = "data/" + string(text) + ".pdb";
         cout << "Input file is currently " << file_in << endl;
      };

   auto out = input_box("Output path");
   out.second->set_text("figures/");
   out.second->on_text =
      [input = out.second.get()](std::string_view text)
      {
         file_out = text;
         cout << "Output file is currently " << file_out << endl;
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
         cout << "Selected " << select << endl;
         if (select == "Radial strategy") {
            setting::grid::psc = setting::grid::RadialStrategy;
         } else if (select == "Axes strategy") {
            setting::grid::psc = setting::grid::AxesStrategy;
         } else if (select == "Jan strategy") {
            setting::grid::psc = setting::grid::JanStrategy;
         }
      },
      {
         "Radial strategy",
         "Axes strategy",
         "Jan strategy"
      }

   ).first; // We return only the first, the menu. the second is a shared pointer to the label.
}

auto image_control(view& _view) {
   auto radio_distance = radio_button("Distance histogram");
   auto radio_intensity = radio_button("Scattering intensity");
   radio_distance.select(true);

   radio_distance.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          cout << "Clicked on second button!" << endl;
          *current_image = *image_distance;
          current_image_selection = 1;
          _view.refresh();
      }
   };

   radio_intensity.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          cout << "Clicked on first button!" << endl;
          *current_image = *image_intensity;
          current_image_selection = 2;
          _view.refresh();
      }
   };

   return group("Change the shown figure",
            margin({10, 10, 20, 20},
               top_margin(25,
                  htile
                  (
                     top_margin(10, align_left(radio_intensity)),
                     top_margin(10, align_left(radio_distance))
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
            Protein protein(file_in);
            protein.generate_new_hydration();
            shared_ptr<ScatteringHistogram> d = protein.get_histogram();

            // Distance plot
            PlotDistance d_plot(d);
            string path = file_out + "distances.png";
            d_plot.save(path); 
            image_distance = share(image("figures/distances.png"));

            // Debye scattering intensity plot
            PlotIntensity i_plot(d);
            path = file_out + "intensity.png";
            i_plot.save(path);
            image_intensity = share(image("figures/intensity.png"));

            if (current_image_selection == 1) {*current_image = *image_distance;} 
            else if (current_image_selection == 2) {*current_image = *image_intensity;} 

            _view.refresh();
         }
      };
   return lbutton;
}

auto make_controls(view& _view) {
   image_intensity = share(image("intensity.png"));
   image_distance = share(image("distances.png"));

   current_image = share(image("distances.png")); //! Must be assigned this way - otherwise compiler optimizations may remove later assignments
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
