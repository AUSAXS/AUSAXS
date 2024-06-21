#pragma once

#include <fitter/FitterFwd.h>

#include <utility/observer_ptr.h>

#include <string>

namespace fitter {
    struct IFit {
        virtual ~IFit() = default;            
        /**
         * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom.
         */ 
        void add_fit(observer_ptr<Fitter> fit) noexcept;

        /**
         * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom. 
         */
        void add_fit(observer_ptr<IFit> fit) noexcept;

        /**
         * @brief Add plots to this fit.
         */
        void add_plots(observer_ptr<Fitter> fitter);

        /**
         * @brief Get a string representation of this object. 
         */
        [[nodiscard]] virtual std::string to_string() const noexcept;
    }
}