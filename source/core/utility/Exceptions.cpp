/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Exceptions.h>
#include <utility/Console.h>

using namespace ausaxs::except;

base::base(const char* msg) : msg(msg) {
    console::print_critical(msg);
}

base::base(const std::string msg) : msg(msg) {
    console::print_critical(msg);
}

const char* base::what() const noexcept {return msg.data();}