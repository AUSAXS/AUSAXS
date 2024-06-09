#pragma once

#include <iostream>
#if defined(_WIN32)
    #include <windows.h>
    #include <map>
#endif

namespace console {
    namespace color {
        typedef int color;
        color black = 0;
        color red = 1;
        color green = 2;
        color yellow = 3;
        color blue = 4;
        color magenta = 5;
        color cyan = 6;
        color lightgray = 7;
        color darkgray = 60;
        color lightred = 61;
        color lightgreen = 62;
        color lightyellow = 63;
        color lightblue = 64;
        color lightmagenta = 65;
        color lightcyan = 66;
        color white = 67;
    }

    #if defined(_WIN32)
        struct detail {
            inline static std::map<color::color, DWORD> foregroundmap = {
                {color::black, 0},
                {color::blue, FOREGROUND_BLUE},
                {color::green, FOREGROUND_GREEN},
                {color::red, FOREGROUND_RED},
                {color::cyan, FOREGROUND_GREEN | FOREGROUND_BLUE},
                {color::magenta, FOREGROUND_RED | FOREGROUND_BLUE},
                {color::lightgray, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE},
                {color::darkgray, FOREGROUND_INTENSITY},
                {color::lightblue, FOREGROUND_BLUE | FOREGROUND_INTENSITY},
                {color::lightgreen, FOREGROUND_GREEN | FOREGROUND_INTENSITY},
                {color::lightred, FOREGROUND_RED | FOREGROUND_INTENSITY},
                {color::lightcyan, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY},
                {color::lightmagenta, FOREGROUND_RED | FOREGROUND_BLUE | FOREGROUND_INTENSITY},
                {color::yellow, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY},
                {color::white, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY}
            };

            inline static std::map<color::color, DWORD> backgroundmap = {
                {color::black, 0},
                {color::blue, BACKGROUND_BLUE},
                {color::green, BACKGROUND_GREEN},
                {color::red, BACKGROUND_RED},
                {color::cyan, BACKGROUND_GREEN | BACKGROUND_BLUE},
                {color::magenta, BACKGROUND_RED | BACKGROUND_BLUE},
                {color::lightgray, BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE},
                {color::darkgray, BACKGROUND_INTENSITY},
                {color::lightblue, BACKGROUND_BLUE | BACKGROUND_INTENSITY},
                {color::lightgreen, BACKGROUND_GREEN | BACKGROUND_INTENSITY},
                {color::lightred, BACKGROUND_RED | BACKGROUND_INTENSITY},
                {color::lightcyan, BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY},
                {color::lightmagenta, BACKGROUND_RED | BACKGROUND_BLUE | BACKGROUND_INTENSITY},
                {color::yellow, BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_INTENSITY},
                {color::white, BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY}
            };
        };
    #endif

    [[maybe_unused]] static void print(std::string message, color::color foreground) {
        #if defined(_WIN32)
            CONSOLE_SCREEN_BUFFER_INFO Info;
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            GetConsoleScreenBufferInfo(hConsole, &Info);
            WORD defaults = Info.wAttributes;
            SetConsoleTextAttribute(hConsole, detail::foregroundmap.at(foreground));
            std::cout << message << std::endl;
            SetConsoleTextAttribute(hConsole, defaults);
        #else
            std::cout << "\033[" << (30 + foreground) << "m"
                        << message
                        << "\033[0m"
                        << std::endl;
        #endif
    }

    [[maybe_unused]] static void print(std::string message, color::color foreground, color::color background) {
        #if defined(_WIN32)
            CONSOLE_SCREEN_BUFFER_INFO Info;
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            GetConsoleScreenBufferInfo(hConsole, &Info);
            WORD defaults = Info.wAttributes;
            SetConsoleTextAttribute(hConsole, detail::foregroundmap.at(foreground) | detail::backgroundmap.at(background));
            std::cout << message << std::endl;
            SetConsoleTextAttribute(hConsole, defaults);
        #else
            std::cout << "\033[" << (30 + foreground) << "m" 
                        << "\033[" << (40 + background) << "m"
                        << message
                        << "\033[0m" 
                        << std::endl;
        #endif
    }
}