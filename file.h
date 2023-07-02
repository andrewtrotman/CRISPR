/*
	FILE.H
	------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief OS agnostic file reading
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include "file_read_only.h"

#include <vector>
#include <iostream>

/*
	CLASS FILE
	----------
*/
class file
	{
	public:
		/*
			FILE::READ_ENTIRE_FILE()
			------------------------
		*/
		/*!
			@brief Read the contents of file filename into the std::string into.
			@details Because into is a string it is naturally '\0' terminated by the C++ std::string class.
			@param filename [in] The path of the file to read.
			@param into [out] The std::string to write into.  This string will be re-sized to the size of the file.
			@return The size of the file in bytes
		*/
		static size_t read_entire_file(const std::string &filename, std::string &into);
	};

	/*
		READ_GUIDES()
		-------------
	*/
	/*!
		@brief Read a load of 20-mers from disk, encoded one per line (everything after the 20th character on each line is ignored)
		@param filename [in] The name of the file to read.
		@param pack_20mer [in] A functor that will pack a read sequence into a 64-bit integer.
		@param map [in] if true mmap the file, if false then use OS file I/O methods.
		@returns A sorted vector of the packed sequences once read from disk.
	*/
	template <typename PACKER>
	std::vector<uint64_t> read_guides(const std::string &filename, PACKER pack_20mer, bool map = true)
		{
		std::vector<uint64_t> packed_guides;
		std::string data;
		file_read_only memory_map;
		char *guide;
		const uint8_t *address_in_memory;
		size_t length;

		if (map)
			{
			memory_map.open(filename);
			length = memory_map.read_entire_file(address_in_memory);
			}
		else
			{
			length = file::read_entire_file(filename, data);
			address_in_memory = (const uint8_t *)&data[0];
			}
		guide = (char *)address_in_memory - 1;

		if (length == 0)
			{
			std::cerr << "Error opening guide file: " << filename << std::endl;
			exit(1);
			}

		do
			{
			guide++;
			packed_guides.push_back(pack_20mer(guide));
			guide = strchr(guide + 20, '\n');
			}
		while (((uint8_t *)guide - (uint8_t *)address_in_memory) < length - 1);

		std::sort(packed_guides.begin(), packed_guides.end());

		return packed_guides;
		}
