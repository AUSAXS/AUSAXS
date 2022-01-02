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
    auto im_intensity = image{"intensity.png"};
    auto im_distance = image{"distances.png"};
    auto radio_intensity = radio_button("first option");
    auto radio_distance = radio_button("second option");
    auto radio_3 = radio_button("third option");
    radio_intensity.select(true);

    auto image_shown =  share(
                            align_center_bottom(
                                fixed_size(
                                    {1000, 700},
                                    im_intensity
                                )
                            )
                        );

    view_.add(image_shown);
    auto image_box =    align_center_bottom(
                            fixed_size(
                                {1000, 700},
                                box(bkd_color_accent)
                            )
                        );

    radio_intensity.on_click = [&view_, &image_shown, &im_intensity] (bool) mutable {
        // if (!radio_distance.value()) {return;}
        cout << "Clicked on first button!" << endl;
        view_.remove(image_shown);
        
        image_shown =   share(
                            align_center_bottom(
                                fixed_size(
                                    {1000, 700},
                                    im_intensity
                                )
                            )
                        );
        view_.add(image_shown);

        // im_shown = im_intensity;
        // view_.refresh(im_shown);
    };

    radio_distance.on_click = [&view_, &image_shown, &im_distance] (bool) mutable {
        // if (!radio_intensity.value()) {return;}
        cout << "Clicked on second button!" << endl;
        view_.remove(image_shown);
    
        image_shown =   share(
                            align_center_bottom(
                                fixed_size(
                                    {1000, 700},
                                    im_distance
                                )
                            )
                        );
        view_.add(image_shown);
        // im_shown = im_distance;
        // view_.refresh(im_shown);
    };

   radio_3.on_click = [&view_, &image_shown] (bool) mutable {
        cout << "Clicked on third button!" << endl;
        view_.remove(image_shown);
   };

   auto radio_buttons =
        group("Radio Buttons",
        margin({10, 10, 20, 20},
            top_margin(25,
                htile(
                    top_margin(10, align_left(radio_intensity)),
                    top_margin(10, align_left(radio_distance)),
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

    view_.content(
        make_controls(view_),
        background
    );

   _app.run();
   return 0;
}