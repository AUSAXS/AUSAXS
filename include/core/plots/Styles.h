#pragma once

#include <string>

namespace ausaxs::style {
    typedef std::string LineStyle;
    typedef std::string DrawStyle;
	typedef std::string Color;

    struct color {
        inline static const Color black =  "k";
        inline static const Color blue =   "tab:blue";
        inline static const Color orange = "tab:orange";
        inline static const Color green =  "tab:green";
        inline static const Color red =    "tab:red";
        inline static const Color purple = "tab:purple";
        inline static const Color brown =  "tab:brown";
        inline static const Color pink =   "tab:pink";
        inline static const Color gray =   "tab:gray";
        inline static const Color olive =  "tab:olive";
        inline static const Color cyan =   "tab:cyan";


        static Color next() {
            static unsigned int i = 0;
            switch (++i % 11) {
                case 1:  return orange;
                case 2:  return blue;
                case 3:  return green;
                case 4:  return red;
                case 5:  return purple;
                case 6:  return brown;
                case 7:  return pink;
                case 8:  return gray;
                case 9:  return olive;
                case 10: return cyan;
                default: return black;
            }
        }        
    };

    namespace color_map {
        struct ColorMap {
            ColorMap(unsigned int n) : n(n) {}
            virtual std::string next() = 0;
            unsigned int i = 0, n;
        };

        struct Rainbow : public ColorMap {
            using ColorMap::ColorMap;
            std::string next() override;
        };

        struct RedBlue : public ColorMap {
            using ColorMap::ColorMap;
            std::string next() override;
        };
    }

    struct line {
        inline static const LineStyle solid = "-";
        inline static const LineStyle dashed = "--";
        inline static const LineStyle dotted = ":";
        inline static const LineStyle dashdot = "-.";
    };

    struct draw {
        inline static const DrawStyle line = "line";
        inline static const DrawStyle hist = "hist";
        inline static const DrawStyle points = "points";
        inline static const DrawStyle errors = "errors";
    };
}