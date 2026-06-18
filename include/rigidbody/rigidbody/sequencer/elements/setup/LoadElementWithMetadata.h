// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/setup/LoadElement.h>

#include <cstdint>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief A LoadElement that additionally retains the per-atom PDB metadata needed to draw a structure preview (residue number + Cα flag).
     */
    class LoadElementWithMetadata : public LoadElement {
        public:
            LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names = {});
            LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& split, const std::vector<std::string>& body_names = {});
            LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names = {});
            ~LoadElementWithMetadata() override;

            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

            // Per-base-atom metadata from the most recent parse, aligned with the molecule's base
            // atoms in body order (before symmetry copies are realized). Stored as statics rather
            // than on the molecule/setup since only one sequencer is ever parsed at a time, and the
            // data is consumed immediately afterwards by the preview API.
            inline static std::vector<int> residue_seq;
            inline static std::vector<std::uint8_t> is_ca;

        private:
            // re-read the resolved files and refresh the static metadata above
            void capture_metadata();
    };
}
