// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/pdb/Record.h>
#include <utility/Utility.h>

#include <string>
#include <vector>

namespace ausaxs::io::pdb {
	class Header : Record {
		public: 
			Header() = default;
			Header(const Header& rhs) = default;
			Header(Header&& rhs) noexcept = default;
			Header &operator=(const Header& rhs) = default;
			Header &operator=(Header&& rhs) noexcept = default;
			~Header() override = default;

			/**
			 * @brief Get the RecordType of this object.
			 */
			RecordType get_type() const override;

			/**
			 * @brief Parse a .pdb format header string. This is equivalent to the add method.
			 * @param s the .pdb format header string.
			 */
			void parse_pdb(const std::string& s) override;

			/**
			 * @brief Get the .pdb format representation of this Header. This is equivalent to the get method.
			 * @return the .pdb format header string. 
			 */
			std::string as_pdb() const override;

			/**
			 * @brief Add a header line to the internal storage of this Header. 
			 * @param s the header line. 
			 */
			void add(const std::string& s);

			/**
			 * @brief Get the .pdb format representation of this Header.
			 * @return the .pdb format header string. 
			 */
			std::string get() const;

			/**
			 * @brief Remove all records of a given type. 
			 */
			void remove(const std::string& type);

			/**
			 * @brief Get the number of header lines.
			 */
			std::size_t size() const;

			bool operator==(const Header& rhs) const;

		private: 
			std::vector<std::string> contents;
	};
}