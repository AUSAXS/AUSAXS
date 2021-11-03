#pragma once
#include <TStyle.h>
#include <TROOT.h>

double labelsize = 0.06;
double titlesize = 0.07;
double xlabeloffset = 0.7;
double ylabeloffset = 0.65;

void setup_style() {
    gStyle->SetLabelSize(labelsize, "X");
    gStyle->SetLabelSize(labelsize, "Y");
    gStyle->SetLabelSize(labelsize, "Z");
    gStyle->SetTitleSize(titlesize, "X");
    gStyle->SetTitleSize(titlesize, "Y");
    gStyle->SetTitleOffset(xlabeloffset, "X");
    gStyle->SetTitleOffset(ylabeloffset, "Y");
    gStyle->SetPadBottomMargin(0.13);
    // gStyle->SetTitleXSize(0.04);
    // gStyle->SetTitleYSize(0.04);
    // gStyle->SetTickLength(0.04);
    gStyle->SetLineStyleString(11, "20 10"); // smaller dashes than the standard
    gStyle->SetLineStyleString(12, "20 20"); // smaller, more spread out dashes than the standard
    gStyle->SetPalette(kViridis); // set the global color scheme of figures
    gStyle->SetOptStat(0); // hide legends
    gStyle->SetOptTitle(0); // hide titles
    gROOT->ForceStyle();
}
