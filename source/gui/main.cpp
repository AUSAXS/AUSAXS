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
auto constexpr bkd_color = rgba(35, 35, 37, 255);
auto background = box(bkd_color);

constexpr auto bred     = colors::red.opacity(0.4);
constexpr auto bgreen   = colors::green.level(0.7).opacity(0.4);
constexpr auto bblue    = colors::blue.opacity(0.4);
constexpr auto brblue   = colors::royal_blue.opacity(0.4);
constexpr auto pgold   = colors::gold.opacity(0.8);

auto make_buttons(view& view_) {
   auto mbutton         = button("Momentary Button");
   auto tbutton         = toggle_button("Toggle Button", 1.0, bred);
   auto lbutton         = share(latching_button("Latching Button", 1.0, bgreen));
   auto reset           = button("Clear Latch", icons::lock_open, 1.0, bblue);
   auto note            = button(icons::cog, "Setup", 1.0, brblue);
   auto prog_bar        = share(progress_bar(rbox(colors::black), rbox(pgold)));
   auto prog_advance    = button("Advance Progress Bar");

   reset.on_click =
      [lbutton, &view_](bool) mutable
      {
         lbutton->value(0);
         view_.refresh(*lbutton);
      };

   prog_advance.on_click =
      [prog_bar, &view_](bool) mutable
      {
         auto val = prog_bar->value();
         if (val > 0.9)
            prog_bar->value(0.0);
         else
            prog_bar->value(val + 0.125);
         view_.refresh(*prog_bar);
      };

   return
      margin({ 20, 0, 20, 20 },
         vtile(
            top_margin(20, mbutton),
            top_margin(20, tbutton),
            top_margin(20, hold(lbutton)),
            top_margin(20, reset),
            top_margin(20, note),
            top_margin(20, vsize(25, hold(prog_bar))),
            top_margin(20, prog_advance)
         )
      );
}

auto make_controls(view& view_) {
   auto im_intensity = image{"intensity.png"};
   auto im_distance = image{"distances.png"};
   auto radio_intensity = radio_button("Scattering histogram");
   auto radio_distance = radio_button("Distance histogram");
   radio_intensity.select(true);

   auto im_shown = im_intensity;
   radio_intensity.on_click = [&radio_intensity, &view_, &im_intensity, &im_shown] (bool) mutable {
      cout << "Clicked on intensity button!" << endl;
      // im_shown = im_intensity;
      // view_.refresh(im_shown);
   };

   radio_distance.on_click = [&view_, &im_distance, &im_shown] (bool) mutable {
      cout << "Clicked on distance button!" << endl;
      // im_shown = im_distance;
      // view_.refresh(im_shown);
   };

   auto  radio_buttons =
         group("Radio Buttons",
            margin({10, 10, 20, 20},
               top_margin(25,
                  vtile(
                     top_margin(10, align_left(radio_intensity)),
                     top_margin(10, align_left(radio_distance))
                  )
               )
            )
         );

   return
      vtile(
         htile(
            margin({ 20, 20, 20, 20 }, radio_buttons)
         ),
         im_shown
      );
}

int main(int argc, char* argv[]) {
   app _app(argc, argv, "SAXS", "com.saxs.gui");
   window _win(_app.name(), 15, rect{20, 20, 1000, 1500});
   _win.on_close = [&_app]() { _app.stop(); };

   view view_(_win);

   view_.content(
      make_controls(view_),
      background
   );

   _app.run();
   return 0;
}