#pragma once

#include <TROOT.h>

/**
 * @brief Defines the plot style used for a single plot.
 */
class PlotOptions {
  public:
    int color = kBlack;
    double alpha = 1;
    int marker_style = 7;
    bool line = true;
    bool markers = false;
    bool use_existing_axes = true;
};