// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/detail/header/data/HeaderData.h>

#include <array>

namespace ausaxs::em::detail::header {
    /**
     * @brief The 1024-byte header of IMOD REC files, as specified in https://bio3d.colorado.edu/imod/doc/mrc_format.txt. 
     *        It assumes the endian of the data is the same as the system. 
     */
    struct RECData {
        RECData();
        ~RECData() = default;

        //! members CANNOT be reordered!
        int nx = 0;                         // Number of columns. Spans the axis given by mapc.
        int ny = 0;                         // Number of rows. Spans the axis given by mapr.
        int nz = 0;                         // Number of sections. Spans the axis given by maps.
        int mode = 0;                       // Bit format of the data.
        int nxstart = 0;                    // Location of first column in the unit cell.
        int nystart = 0;                    // Location of first row in the unit cell.
        int nzstart = 0;                    // Location of first section in the unit cell. 
        int mx = 0;                         // Sampling rate along the x-axis in the unit cell.
        int my = 0;                         // Sampling rate along the y-axis in the unit cell.
        int mz = 0;                         // Sampling rate along the z-axis in the unit cell.
        float cella_x = 0;                  // Unit cell size in Ångström along x-axis. The map may store only part of it.
        float cella_y = 0;                  // Unit cell size in Ångström along y-axis. The map may store only part of it.
        float cella_z = 0;                  // Unit cell size in Ångström along z-axis. The map may store only part of it.
        float cellb_alpha = 0;              // Cell angle alpha in degrees.
        float cellb_beta = 0;               // Cell angle beta in degrees.
        float cellb_gamma = 0;              // Cell angle gamma in degrees.
        int mapc = 0;                       // Axis corresponding to columns (1,2,3 for X,Y,Z).
        int mapr = 0;                       // Axis corresponding to rows (1,2,3 for X,Y,Z).
        int maps = 0;                       // Axis corresponding to sections (1,2,3 for X,Y,Z).
        float dmin = 0;                     // Minimum pixel value.
        float dmax = 0;                     // Maximum pixel value.
        float dmean = 0;                    // Mean pixel value.
        int ispg = 0;                       // Space group number.
        int nsymbt = 0;                     // Size of extended header in bytes.
        short int creatid = 0;              // Creator ID.
        std::array<char, 6> extra1{};
        std::array<char, 4> exttyp{};       // Type of extended header.
        int nversion = 0;                   // Version number of the format.
        std::array<char, 16> extra2{};
        short int nint = 0;                 // Number of integers per section.
        short int nreal = 0;                // Number of reals per section.
        std::array<char, 20> extra3{};
        int imodstamp = 0;                  // IMOD stamp.
        int imodflags = 0;                  // IMOD flags.
        short int idtype = 0;               // ID type.
        short int lens = 0;                 // Lens type.
        short int nd1 = 0;                  // ND1.
        short int nd2 = 0;                  // ND2.
        short int vd1 = 0;                  // VD1.
        short int vd2 = 0;                  // VD2.
        std::array<float, 6> tiltangles{};  // Tilt angles.
        float origin_x = 0;
        float origin_y = 0;
        float origin_z = 0;
        std::array<char, 4> map{};
        std::array<char, 4> stamp{};
        float rms = 0;
        int nlabl = 0;
        std::array<char, 800> label{};
    };
    static_assert(sizeof(RECData) == 1024, "em::detail::header::RECData: Size of RECData is wrong.");
}