/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/Flags.h>

using namespace ausaxs::settings;

bool flags::data_rebin = false;
char flags::last_parsed_unit = ' ';
bool flags::init_histogram_manager = true;