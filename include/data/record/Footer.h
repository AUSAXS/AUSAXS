#pragma once

#include <data/record/Record.h>
#include <utility/Utility.h>

#include <string>
#include <vector>

namespace data::record {
	class Footer : Record {
		public: 
			Footer() noexcept;

			~Footer() override;

			/**
			 * @brief Get the RecordType of this object.
			 */
			RecordType get_type() const override;

			/**
			 * @brief Parse a .pdb format Footer string. This is equivalent to the add method.
			 * @param s the .pdb format Footer string.
			 */
			void parse_pdb(const std::string& s) override;

			/**
			 * @brief Get the .pdb format representation of this Footer. This is equivalent to the get method.
			 * @return the .pdb format Footer string. 
			 */
			std::string as_pdb() const override;

			/**
			 * @brief Add a Footer line to the internal storage of this Footer. 
			 * @param s the Footer line. 
			 */
			void add(const std::string& s);

			/**
			 * @brief Remove all records of a given type. 
			 */
			void remove(const std::string& type);

			/**
			 * @brief Get the .pdb format representation of this Footer.
			 * @return the .pdb format Footer string. 
			 */
			std::string get() const;

			/**
			 * @brief Get the number of header lines.
			 */
			unsigned int size() const;

			bool operator==(const Footer& rhs) const;

		private: 
			std::vector<std::string> contents;
	};
}