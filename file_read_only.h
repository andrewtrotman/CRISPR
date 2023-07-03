/*
	FILE_READ_ONLY.H
	----------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief OS agnostic mmap() stuff for reading files
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>

#ifdef _MSC_VER
	#include <Windows.h>
#else
	#include <unistd.h>
	#include <sys/mman.h>
#endif

#include <string>

/*
	CLASS FILE_READ_ONLY
	--------------------
*/
/*!
	@brief A read_only file object, the memory was probably allocated with mmap() and needs deallocating accordingly
*/
class file_read_only
	{
	private:
#ifdef _MSC_VER
		HANDLE hFile;							///< The file being mapped
		HANDLE hMapFile;						///< The mapping of that file
#endif
		const void *file_contents;				///< The contents of the file.
		size_t size;							///< The size of the file.

	public:
		/*
			FILE::FILE_READ_ONLY::FILE_READ_ONLY()
			--------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		file_read_only():
			file_contents(nullptr),
			size(0)
			{
			}

		/*
			FILE_READ_ONLY::~FILE_READ_ONLY()
			---------------------------------
		*/
		/*!
			@brief Destructor
		*/
		~file_read_only()
			{
			#ifdef _MSC_VER
				UnmapViewOfFile((void *)file_contents);
				CloseHandle(hMapFile); // close the file mapping object
				CloseHandle(hFile);   // close the file itself
			#else
				munmap((void *)file_contents, size);
			#endif
			}

		/*
			FILE_READ_ONLY::OPEN()
			----------------------
		*/
		/*!
			@brief Open and read the file into memory
			@param filename [in] The name of the file to read
			@return The size of the file
		*/
		size_t open(const std::string &filename);

		/*
			FILE_READ_ONLY::READ_ENTIRE_FILE()
			----------------------------------
		*/
		/*!
			@brief Return the contents and length of the file.
			@param into [out] The pointer to write into.
			@return The size of the file in bytes
		*/
		size_t read_entire_file(const uint8_t *&into) const
			{
			into = reinterpret_cast<const uint8_t *>(file_contents);
			return size;
			}
	};
