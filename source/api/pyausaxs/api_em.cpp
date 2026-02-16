// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_pyausaxs.h>
#include <api/ObjectStorage.h>
#include <io/Reader.h>
#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/detail/SimpleExvModel.h>
#include <fitter/SmartFitter.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;

// #include <em/ImageStack.h>
// int map_read(
//     const char* filename,
//     int* status
// ) {return execute_with_catch([&]() {
//     return map_id;
// }, status);}

// void map_get_slice(
//     int map_id,
//     double z_position,
//     double** slice_data,
//     int* width, int* height,
//     int* status
// );

// int map_fit(
//     int map_id, int data_id,
//     int* status
// );