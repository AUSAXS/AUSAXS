// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#pragma once

#include <string_view>

namespace ausaxs::logging {
    /**
     * @brief Write a message to the log.
     */
    void log(std::string_view msg);

    /**
     * @brief Replicate a console message to the log.
     *        This will prepend the message with "CONSOLE: ".
     */
    void log_console(std::string_view msg);

    /**
     * @brief Write a critical message to the log.
     *        This will immediately flush all pending messages to the log.
     */
    void log_critical(std::string_view msg);

    /**
     * @brief Start a new log named `name`_log.
     */
    void start(std::string_view name);
}