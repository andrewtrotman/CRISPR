/*
	FILE_READ_ONLY.CPP
	------------------
*/
#include "file_read_only.h"

/*
	FILE_READ_ONLY::OPEN()
	----------------------
*/
size_t file_read_only::open(const std::string &filename)
	{
	#ifdef _MSC_VER
		hFile = CreateFile(filename.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_READONLY, NULL);
		if (hFile == INVALID_HANDLE_VALUE)
			return 0;

		hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
		if (hMapFile == NULL)
			{
			CloseHandle(hFile);
			return 0;
			}

		void *lpMapAddress = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);
		if (lpMapAddress == NULL)
			{
			CloseHandle(hFile);
			CloseHandle(hMapFile);
			return 0;
			}

		file_contents = (uint8_t *)lpMapAddress;

		DWORD high;
		DWORD low = GetFileSize(hFile, &high);

		size = ((uint64_t)high << (uint64_t)32) + (uint64_t)low;

		return size;
	#else
		/*
			Open the file
		*/
		int reader;

		if ((reader = ::open(filename.c_str(), O_RDONLY)) < 0)
			return 0;

		/*
			Find out how large it is
		*/
		struct stat statistics;
		if (fstat(reader, &statistics) != 0)
			{
			close(reader);
			return 0;
			}

		/*
			Allocate space for it and load it
		*/
		#ifdef __APPLE__
			file_contents = (uint8_t *)mmap(nullptr, statistics.st_size, PROT_READ, MAP_PRIVATE, reader, 0);
		#else
			file_contents = (uint8_t *)mmap(nullptr, statistics.st_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, reader, 0);
		#endif

		/*
			Close the file
		*/
		close(reader);

		if (file_contents == nullptr)
			return 0;

		/*
			Remember the file size
		*/
		size = statistics.st_size;

		return size;
	#endif
	}
