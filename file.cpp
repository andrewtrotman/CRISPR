/*
	FILE.CPP
	--------
	Copyright (c) 2023 Andrew Trotman
*/
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "file.h"

/*
	FILE::READ_ENTIRE_FILE()
	------------------------
	This uses a combination of "C" FILE I/O and C++ strings in order to copy the contents of a file into an internal buffer.
	There are many different ways to do this, but this is the fastest according to this link: http://insanecoding.blogspot.co.nz/2011/11/how-to-read-in-file-in-c.html
	Note that there does not appear to be a way in C++ to avoid the initialisation of the string buffer.

	Returns the length of the file in bytes - which is also the size of the string buffer once read.
*/
size_t file::read_entire_file(const std::string &filename, std::string &into)
	{
	FILE *fp;
	// "C" pointer to the file
#ifdef _MSC_VER
	struct __stat64 details;				// file system's details of the file
#else
	struct stat details;				// file system's details of the file
#endif
	size_t file_length = 0;			// length of the file in bytes

	/*
		Fopen() the file then fstat() it.  The alternative is to stat() then fopen() - but that is wrong because the file might change between the two calls.
	*/
	if ((fp = fopen(filename.c_str(), "rb")) != nullptr)
		{
#ifdef _MSC_VER
		if (_fstat64(fileno(fp), &details) == 0)
#else
		if (fstat(fileno(fp), &details) == 0)
#endif
			if ((file_length = details.st_size) != 0)
				{
				into.resize(file_length);
				if (fread(&into[0], details.st_size, 1, fp) != 1)
					into.resize(0);				// LCOV_EXCL_LINE	// happens when reading the file_size buyes failes (i.e. disk or file failure).
				}
		fclose(fp);
		}

	return file_length;
	}
