#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/ExvTable.h>
#include <form_factor/ExvFormFactor.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::form_factor {
    class ExvTableManager {
        public:
            static observer_ptr<const constants::exv::detail::ExvSet> get_current_exv_table();
            static const detail::ExvFormFactorSet& get_current_exv_form_factor_set();
            void set_custom_exv_table(const constants::exv::detail::ExvSet& set);

        private:
            bool _use_custom_exv_table = false;
            static std::unique_ptr<constants::exv::detail::ExvSet> custom_exv_tables;
    };
}