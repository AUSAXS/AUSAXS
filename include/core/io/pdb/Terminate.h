// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/pdb/Record.h>

namespace ausaxs::io::pdb {
	/**
	 * @brief A representation of a Terminate record. 
	 */
	class Terminate : public Record {
		public: 
			Terminate() = default;
			Terminate(const Terminate& rhs) = default;
			Terminate(Terminate&& rhs) noexcept = default;
			Terminate &operator=(const Terminate& rhs) = default;
			Terminate &operator=(Terminate&& rhs) noexcept = default;
			~Terminate() override = default;

			/**
			 * @brief Constructor.
			 * 
			 * @param serial Serial number of the record. 
			 * @param resName Residue name. 
			 * @param chainID The chain identifier. 
			 * @param resSeq The residue sequence identifier. 
			 * @param iCode iCode. 
			 */
			Terminate(int serial, const std::string& resName, char chainID, int resSeq, const std::string& iCode);

			/**
			 * @brief Get the RecordType of this object.
			 */
			RecordType get_type() const override;

			/**
			 * @brief Parse a .pdb format terminate string.
			 * 
			 * @param s the .pdb format terminate string.
			 */
			void parse_pdb(const std::string& s) override;

			/**
			 * @brief Get the .pdb format representation of this Terminate. This is equivalent to the get method.
			 * @return the .pdb format header string. 
			 */
			std::string as_pdb() const override;

			/**
			 * @brief Get the .pdb format representation of this Terminate.
			 */
			std::string get() const;
			
			/**
			 * @brief Set the serial of this record.
			 */
			void set_serial(int serial);

			bool operator==(const Terminate& rhs) const;

			int serial, resSeq;
			std::string resName, iCode;
			char chainID;
	};
}