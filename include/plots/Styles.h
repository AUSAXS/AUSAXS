#pragma once

#include <string>

struct style {
    typedef std::string LineStyle;
    typedef std::string DrawStyle;
	typedef std::string Color;

    struct color {
        inline static Color black = "k";
        inline static Color blue = "tab:blue";
        inline static Color orange = "tab:orange";
        inline static Color green = "tab:green";
        inline static Color red = "tab:red";
        inline static Color purple = "tab:purple";
        inline static Color brown = "tab:brown";
        inline static Color pink = "tab:pink";
        inline static Color gray = "tab:gray";
        inline static Color olive = "tab:olive";
        inline static Color cyan = "tab:cyan";
    };

    struct line {
        inline static LineStyle solid = "-";
        inline static LineStyle dashed = "--";
        inline static LineStyle dotted = ":";
        inline static LineStyle dashdot = "-.";
    };

    struct draw {
        inline static DrawStyle line = "line";
        inline static DrawStyle hist = "hist";
        inline static DrawStyle points = "points";
        inline static DrawStyle errors = "errors";
    };
};