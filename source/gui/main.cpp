/*=============================================================================
   Copyright (c) 2016-2020 Joel de Guzman

   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/
#include <elements.hpp>
#include <algorithm>
#include <random>
#include <iostream>

using namespace cycfi::elements;
using std::cout, std::endl; 

// Main window background color
auto constexpr bkd_color_accent = rgba(55, 55, 57, 255);
auto constexpr bkd_color = rgba(35, 35, 37, 255);
auto background = box(bkd_color);

// The different image types which can be plotted
std::shared_ptr<cycfi::elements::image> image_intensity, image_distance; 

// The currently shown image
std::shared_ptr<cycfi::elements::image> current_image;

auto image_control(view& _view) {
   auto radio_intensity = radio_button("Scattering intensity");
   auto radio_distance = radio_button("Distance histogram");
   radio_intensity.select(true);

   radio_intensity.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          cout << "Clicked on first button!" << endl;
          *current_image = *image_intensity;
          _view.refresh();
      }
   };

   radio_distance.on_click = [&_view] (bool pressed) mutable {
      if (pressed) {
          cout << "Clicked on second button!" << endl;
          *current_image = *image_distance;
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

auto make_controls(view& _view) {
   image_intensity = share(image("intensity.png"));
   image_distance = share(image("distances.png"));

   current_image = share(image("intensity.png")); //! Must be assigned this way - otherwise compiler optimizations may remove later assignments
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
               htile
               (
                  margin({20, 20, 20, 20}, image_control(_view))
               ),
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
      // make_controls(_view),
      make_controls(_view),
      background
   );

   _app.run();
   return 0;
}
