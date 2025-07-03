// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/detail/header/data/HeaderData.h>

#include <array>

namespace ausaxs::em::detail::header {
    /**
     * @brief The 1024-byte header of MRC files as specified in https://www.ccpem.ac.uk/mrc_format/mrc2014.php.
     *        It assumes the endian of the data is the same as the system. 
     */
    struct MRCData {
        MRCData();
        ~MRCData();

        //! members CANNOT be reordered!
        int nx;             // Number of points along x-axis.
        int ny;             // Number of points along y-axis.
        int nz;             // Number of points along z-axis.
        int mode;           // Bit format of the data.
        int nxstart;        // Location of first column in the unit cell.
        int nystart;        // Location of first row in the unit cell.
        int nzstart;        // Location of first section in the unit cell. 
        int mx;             // Sampling rate along the x-axis in the unit cell.
        int my;             // Sampling rate along the y-axis in the unit cell.
        int mz;             // Sampling rate along the z-axis in the unit cell.
        float cella_x;      // Map size in Ångström along x-axis. 
        float cella_y;      // Map size in Ångström along y-axis.
        float cella_z;      // Map size in Ångström along z-axis.
        float cellb_alpha;  // Cell angle alpha in degrees.
        float cellb_beta;   // Cell angle beta in degrees.
        float cellb_gamma;  // Cell angle gamma in degrees.
        int mapc;           // Axis corresponding to columns (1,2,3 for X,Y,Z).
        int mapr;           // Axis corresponding to rows (1,2,3 for X,Y,Z).
        int maps;           // Axis corresponding to sections (1,2,3 for X,Y,Z).
        float dmin;         // Minimum pixel value.
        float dmax;         // Maximum pixel value.
        float dmean;        // Mean pixel value.
        int ispg;           // Space group number.
        int nsymbt;         // Size of extended header in bytes.
        std::array<char, 8> extra1;
        std::array<char, 4> exttyp;
        int nversion;       // Version number of the format.
        std::array<char, 84> extra2;
        float origin_x;
        float origin_y;
        float origin_z;
        std::array<char, 4> map;
        unsigned int machst;
        float rms;
        int nlabl;
        std::array<char, 800> label;
    };
    static_assert(sizeof(MRCData) == 1024, "em::detail::header::MRCData: Size of MRCData is wrong.");
}