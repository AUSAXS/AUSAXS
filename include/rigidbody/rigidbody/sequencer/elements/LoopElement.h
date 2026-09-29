// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <fitter/FitterFwd.h>
#include <io/IOFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>

#include <atomic>
#include <memory>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief A loop element is a sequence element that repeats whatever is inside it a number of times.
     */
    class LoopElement : public GenericElement {
        friend class OptimizeStepElement;
        public:
            LoopElement(observer_ptr<LoopElement> owner, int repeats);
            ~LoopElement() override;

            virtual std::shared_ptr<fitter::FitResult> execute();

            /**
             * @brief Run an iteration of this loop. 
             */
            void run() override;

            virtual observer_ptr<Rigidbody> _get_rigidbody() const;

            virtual observer_ptr<data::Molecule> _get_molecule() const;

            virtual observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> _get_best_conf() const;
            virtual observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> _get_current_conf() const;

            virtual observer_ptr<const Sequencer> _get_sequencer() const;
            virtual observer_ptr<Sequencer> _get_sequencer();

            std::vector<std::unique_ptr<GenericElement>>& _get_elements();
            int _get_loop_iterations() const;

            observer_ptr<LoopElement> _get_owner() const;

            /**
             * @brief Request that the currently running sequence stops as soon as possible.
             *
             * Thread-safe, and intended to be called from outside the running thread. The flag is checked at the
             * start of every loop iteration, so the iteration in progress is always allowed to finish. The stopped
             * run still completes normally, i.e. the best conformation found so far is restored and fitted.
             * The flag is cleared by Sequencer::execute, so a request made while nothing is running is discarded.
             */
            static void _request_stop();
            static bool _stop_requested();
            static void _clear_stop_request();

            static int _get_current_iteration();
            static int _get_total_iterations();
            /**
             * @brief Recalculate the total number of optimization steps the given element tree will perform.
             *
             * This must be done after the tree is fully built, since a loop does not know its own contents
             * while it is being constructed.
             */
            static void _recount_total_iterations(observer_ptr<LoopElement> root);
            static void _reset_counters();
            static void _reset_named_loops();

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        protected: 
            int iterations = 1;
            std::vector<std::unique_ptr<GenericElement>> elements;

        private:
            observer_ptr<LoopElement> owner;
            inline static int total_loop_count = 0;
            inline static int global_counter = 0;
            inline static std::atomic<bool> stop_flag = false;
    };
}