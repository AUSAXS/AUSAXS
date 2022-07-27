#pragma once

// includes
#include <string>
#include <vector>
#include <data/Record.h>
#include <utility/Utility.h>

class Header : Record {
public: 
	Header() {}

	Header(const Header& header) noexcept : contents(header.contents) {}

	Header(const Header&& header) noexcept : contents(std::move(header.contents)) {}

	~Header() override {}

	/**
	 * @brief Get the RecordType of this object.
	 * @return Record::HEADER
	 */
	RecordType get_type() const override {return HEADER;}

	/**
	 * @brief Parse a .pdb format header string. This is equivalent to the add method.
	 * @param s the .pdb format header string.
	 */
	void parse_pdb(std::string s) override;

	/**
	 * @brief Get the .pdb format representation of this Header. This is equivalent to the get method.
	 * @return the .pdb format header string. 
	 */
	std::string as_pdb() const override;

	/**
	 * @brief Add a header line to the internal storage of this Header. 
	 * @param s the header line. 
	 */
	void add(const std::string s);

	/**
	 * @brief Get the .pdb format representation of this Header.
	 * @return the .pdb format header string. 
	 */
	std::string get() const;

	/**
	 * @brief Remove all records of a given type. 
	 */
	void remove(std::string type);

	/**
	 * @brief Get the number of header lines.
	 */
	size_t size() const;

	Header& operator=(const Header& header);

private: 
	std::vector<std::string> contents;
};