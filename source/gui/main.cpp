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

auto make_controls(view& view_) {
   auto image1 = image{"intensity.png"};
   auto image2 = image{"distances.png"};
   auto radio_1 = radio_button("button 1");
   auto radio_2 = radio_button("button 2");
   auto radio_3 = radio_button("button 3");
   radio_1.select(true);

   auto image_shown = image1;
   auto image_box = (layer(
                            margin(
                              {10, 10, 10, 10},
                              image_shown
                            //   link(image_shown) doesnt build
                           ),
                           box(bkd_color_accent)
                        )
                    );

   radio_1.on_click = [&image_shown, &image1] (bool) mutable {
      cout << "Clicked on first button!" << endl;
    //   image_shown = image1;
   };

   radio_2.on_click = [&image_shown, &image2] (bool) mutable {
      cout << "Clicked on second button!" << endl;
    //   image_shown = image2;
   };
   radio_3.on_click = [&image_shown, &image2] (bool) mutable {
      cout << "Clicked on third button!" << endl;
    //   image_shown = image2;
   };

   auto radio_buttons =
         group("Radio Buttons",
            margin({10, 10, 20, 20},
               top_margin(25,
                  htile(
                     top_margin(10, align_left(radio_1)),
                     top_margin(10, align_left(radio_2)),
                     top_margin(10, align_left(radio_3))
                  )
               )
            )
         );

   return
      vtile(
         htile(
            margin({20, 20, 20, 20}, radio_buttons)
         ),
         image_box
      );
}

int main(int argc, char* argv[]) {
   app _app(argc, argv, "SAXS", "com.saxs.gui");
   window _win(_app.name(), 15, rect{20, 20, 1000, 1500});
   _win.on_close = [&_app]() { _app.stop(); };

   view view_(_win);


   auto image1 = image{"intensity.png"};
   auto image2 = image{"distances.png"};
   auto radio_1 = radio_button("button 1");
   auto radio_2 = radio_button("button 2");
   auto radio_3 = radio_button("button 3");
   radio_1.select(true);

   auto image_shown = image1;
   auto image_box = (layer(
                            margin(
                              {10, 10, 10, 10},
                            //   image_shown
                              link(image_shown) //doesnt build
                           ),
                           box(bkd_color_accent)
                        )
                    );

   radio_1.on_click = [&view_, &image_shown, &image1] (bool release) mutable {
      if (release) {
          cout << "Clicked on first button!" << endl;
          image_shown = image1;
          view_.refresh();
      }
   };

   radio_2.on_click = [&view_, &image_shown, &image2] (bool release) mutable {
      if (release) {
          cout << "Clicked on second button!" << endl;
          image_shown = image2;
          view_.refresh();
      }
   };
   radio_3.on_click = [&view_, &image_shown, &image2] (bool release) mutable {
      if (release) {
          cout << "Clicked on third button!" << endl;
          image_shown = image2;
          view_.refresh();
      }
   };

   auto radio_buttons =
         group("Radio Buttons",
            margin({10, 10, 20, 20},
               top_margin(25,
                  htile(
                     top_margin(10, align_left(radio_1)),
                     top_margin(10, align_left(radio_2)),
                     top_margin(10, align_left(radio_3))
                  )
               )
            )
         );

    auto image = vtile(
         htile(
            margin({20, 20, 20, 20}, radio_buttons)
         ),
         image_box
      );
   view_.content(
      image,
      background
   );

   _app.run();
   return 0;
}
