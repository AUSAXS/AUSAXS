#pragma once

#include <em/detail/header/data/HeaderData.h>

#include <array>

namespace em::detail::header {
    /**
     * @brief The 1024-byte header of IMOD REC files, as specified in https://bio3d.colorado.edu/imod/doc/mrc_format.txt. 
     *        It assumes the endian of the data is the same as the system. 
     */
    struct RECData : public HeaderData {
        RECData() = default;
        ~RECData() = default;

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
        short int creatid;  // Creator ID.
        std::array<char, 6> extra1;
        std::array<char, 4> exttyp; // Type of extended header.
        int nversion;       // Version number of the format.
        std::array<char, 16> extra2;
        short int nint;     // Number of integers per section.
        short int nreal;    // Number of reals per section.
        std::array<char, 20> extra3;
        int imodstamp;      // IMOD stamp.
        int imodflags;      // IMOD flags.
        short int idtype;   // ID type.
        short int lens;     // Lens type.
        short int nd1;      // ND1.
        short int nd2;      // ND2.
        short int vd1;      // VD1.
        short int vd2;      // VD2.
        std::array<float, 6> tiltangles;   // Tilt angles.
        float origin_x;
        float origin_y;
        float origin_z;
        std::array<char, 4> map;
        std::array<char, 4> stamp;
        float rms;
        int nlabl;
        std::array<char, 800> label;
    };
    static_assert(sizeof(RECData) == 1024, "em::detail::header::RECData: Size of RECData is wrong.");
}