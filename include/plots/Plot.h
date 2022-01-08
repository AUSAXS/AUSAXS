#pragma once

#include <string.h>
#include <memory.h>

#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>

class Plot {
    public: 
        Plot() {if (!stylized) stylize();}
        virtual ~Plot() = default;

        virtual void save(const std::string& folder) const = 0;    
    private: 
        inline static bool stylized = false;

        static void stylize(const EColorPalette palette = kViridis) {
            // static double labelsize = 0.06;
            // static double titlesize = 0.07;
            // static double xlabeloffset = 0.7;
            // static double ylabeloffset = 0.65;

            // gStyle->SetLabelSize(labelsize, "X");
            // gStyle->SetLabelSize(labelsize, "Y");
            // gStyle->SetLabelSize(labelsize, "Z");
            // gStyle->SetTitleSize(titlesize, "X");
            // gStyle->SetTitleSize(titlesize, "Y");
            // gStyle->SetTitleOffset(xlabeloffset, "X");
            // gStyle->SetTitleOffset(ylabeloffset, "Y");
            // gStyle->SetPadBottomMargin(0.13);
            // gStyle->SetTitleXSize(0.04);
            // gStyle->SetTitleYSize(0.04);
            // gStyle->SetTickLength(0.04);
            gStyle->SetLineStyleString(11, "20 10"); // smaller dashes than the standard
            gStyle->SetLineStyleString(12, "20 20"); // smaller, more spread out dashes than the standard
            gStyle->SetPalette(palette); // set the global color scheme of figures
            gStyle->SetOptStat(0); // hide legends
            gStyle->SetOptTitle(0); // hide titles
            gROOT->ForceStyle();
            stylized = true;
        }
};