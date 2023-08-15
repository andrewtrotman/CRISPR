/*
	MAIN.CPP
	--------
	Copyright (c) 2023 Andrew Trotman
*/
#include <limits>
#include <vector>
#include <random>
#include <thread>
#include <iostream>
#include <algorithm>

#include "file.h"
#include "finder.h"
#include "hamming_distance.h"
#include "encode_kmer_2bit.h"
#include "encode_kmer_3bit.h"
#include "fast_binary_search.h"

/*
	SELECT_PSEUDO_RANDOM_VECTORS()
	------------------------------
*/
void select_pseudo_random_vectors(std::vector<uint64_t> &selected, const std::vector<uint64_t> &guides)
	{
	uint64_t seed = 17;
	uint64_t m = guides.size() - 1;		// Modulus parameter
	uint64_t a = 3;						// Multiplier term
	uint64_t c = 3;						// Increment term
	uint64_t next = seed;

	for (uint64_t which = 0; which < selected.size(); which++)
		{
		selected[which] = guides[next];
		next = ((next * a) + c) % m;
		}
	}

/*
	SELECT_RANDOM_VECTORS()
	-----------------------
*/
void select_random_vectors(std::vector<uint64_t> &selected, const std::vector<uint64_t> &guides)
	{
	constexpr int seed = 17;
	std::mt19937 random(seed);
	std::uniform_int_distribution<int> distribution(0, guides.size() - 1);

	for (size_t which = 0; which < selected.size(); which++)
		selected[which] = guides[distribution(random)];
	}

/*
	UNIQ()
	------
*/
/*!
	@brief Unique (inplace) the data in a vector and return a new vector of frquencies of each unique member
	@param data [in/out] The data to unique
	@param frequencies[out] The frequencies of each unique member.
*/
void uniq(std::vector<uint64_t> &data, std::vector<uint16_t> &frequencies)
	{
	uint64_t was = 0;
	uint64_t count = 1;

	/*
		Allocate space for the frequencies
	*/
	frequencies.reserve(data.size());

	/*
		Unique the vector
	*/
	size_t index = 1;
	while (index < data.size())
		{
		if (data[index] == data[was])
			count++;
		else
			{
			was++;
			data[was] = data[index];
			frequencies.push_back(count > std::numeric_limits<uint16_t>::max() ? std::numeric_limits<uint16_t>::max() : count);
			count = 1;
			}
		index++;
		}
	/*
		Don't forget the last one
	*/
	was++;
	data[was] = data[index];
	frequencies.push_back(count);

	/*
		Resize the data to just the size of the unique list
	*/
	data.resize(was);
	}

/*
	Allocate space for the final set of results
*/
std::vector<std::vector<sequence_score_pair>> all_positions;

/*
	Command line parameters
*/
enum { FAST_BINARY_SEARCH, FAST_BINARY_SEARCH_AVX512, HAMMING_DISTANCE };
const std::string mode_name[] = { "Binary Search", "Binary Search using AVX512", "Hamming Distance" };
size_t TESTSIZE;
int mode = FAST_BINARY_SEARCH;
std::string guides_filename = "OryzaSativaGuides.txt";
const char *output_filename = nullptr;
size_t thread_count = std::thread::hardware_concurrency();



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
std::vector<uint64_t> read_guides(std::vector<uint16_t> &frequencies, const std::string &filename, PACKER pack_20mer, bool map = true)
	{
	std::vector<uint64_t> packed_guides;
	std::string data;
	JASS::file::file_read_only memory_map;
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
		length = JASS::file::read_entire_file(filename, data);
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

	/*
		Sort and uniq the list
	*/
	std::sort(packed_guides.begin(), packed_guides.end());
	uniq(packed_guides, frequencies);

	return packed_guides;
	}

/*
	USAGE()
	-------
*/
int usage(const char *exename)
	{
	std::cout << "Usage:" << exename << "[-b | -B | -h] [-f<filename>] [-o<filename>] [-t<threadcount>\n";
	std::cout << "       -b for binary search [default]\n";
//	std::cout << "       -B for binary search using AVX512 instructions\n";
	std::cout << "       -h for hamming distance\n";
	std::cout << "       -t<threads> the number of threads to use when searching [default = corecount (inc hyperthreads)]\n";
	std::cout << "       -f<filename> use <filename> as the genome\n";
	std::cout << "       -o<filename> use <filename> as the output file\n";

	return 0;
	}

/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	job workload;
	encode_kmer_2bit packer_2bit;
	encode_kmer_3bit packer_3bit;

	auto start_main = std::chrono::steady_clock::now(); // Start timing

	if (argc <= 4)
		{
		for (int arg = 1; arg < argc; arg++)
			{
			if (std::string(argv[arg]) == "-b")
				mode = FAST_BINARY_SEARCH;
			else if (std::string(argv[arg]) == "-B")
				mode = FAST_BINARY_SEARCH_AVX512;
			else if (std::string(argv[arg]) == "-h")
				mode = HAMMING_DISTANCE;
			else if (strncmp(argv[arg], "-f", 2) == 0)
				{
				if (strlen(argv[arg]) == 2)
					guides_filename = std::string(argv[++arg]);
				else
					guides_filename = std::string(argv[arg] + 2);
				}
			else if (strncmp(argv[arg], "-o", 2) == 0)
				{
				if (strlen(argv[arg]) == 2)
					output_filename = argv[++arg];
				else
					output_filename = argv[arg] + 2;
				}
			else if (strncmp(argv[arg], "-t", 2) == 0)
				{
				if (strlen(argv[arg]) == 2)
					thread_count = atol(argv[++arg]);
				else
					thread_count = atol(argv[arg] + 2);
				}
			else
				{
				std::cout << "unknown parameter:" << argv[arg] << "\n";
				exit(usage(argv[0]));
				}
			}
		}
	else
		exit(usage(argv[0]));

	std::cout << "Processing using: " << mode_name[mode] << "\n";
	std::cout << "Using file      : " << guides_filename << "\n";
	std::cout << "Output file     : " << (output_filename == nullptr ? "<stdout>" : output_filename) << "\n";

	if  (output_filename == nullptr)
		workload.output_file = stdout;
	else if ((workload.output_file = fopen(output_filename, "wb")) == nullptr)
		{
		std::cout << "Cannot open output file:" << output_filename << "\n";
		exit(1);
		}

	/*
		Read the guides from the file.
		NOTE:read_guides() sorts the list after loading and before returning it.
	*/
	auto time_io_start = std::chrono::steady_clock::now(); // Start timing
	if (mode == FAST_BINARY_SEARCH || mode == FAST_BINARY_SEARCH_AVX512)
		workload.genome_guides = read_guides(workload.genome_guide_frequencies, guides_filename, packer_2bit.pack_20mer);
	else
		workload.genome_guides = read_guides(workload.genome_guide_frequencies, guides_filename, packer_3bit.pack_20mer);
	auto time_io_end = std::chrono::steady_clock::now(); // Stop timing
	auto time_io_duration = std::chrono::duration_cast<std::chrono::microseconds>(time_io_end - time_io_start).count();
	std::cout << "Load encode time: " << time_io_duration / 1000000.001 << " seconds" << '\n';

	/*
		Generate some samples (sometimes random, and sometimes not, so keep both versions)
	*/
	if (true)
		{
		TESTSIZE = 10000;
		workload.guide.resize(TESTSIZE);
		workload.guide_frequencies.resize(TESTSIZE);
//		select_random_vectors(workload.guide, workload.genome_guides);
		select_pseudo_random_vectors(workload.guide, workload.genome_guides);
		sort(workload.guide.begin(), workload.guide.end());
		for (size_t which = 0; which < TESTSIZE; which++)
			workload.guide_frequencies[which] = 1;
		}
	else
		{
		TESTSIZE = workload.genome_guides.size();
		workload.guide.resize(TESTSIZE);
		workload.guide_frequencies.resize(TESTSIZE);
		for (size_t which = 0; which < TESTSIZE; which++)
			{
			workload.guide[which] = workload.genome_guides[which];
			workload.guide_frequencies[which] = workload.genome_guide_frequencies[which];
			}
		}

	std::cout << "Loaded " << workload.guide.size() << " test guides, " <<  workload.genome_guides.size() << " genome guides" << '\n';

	/*
		Generate the index over the guides.
	*/
	finder *searcher = nullptr;
	if (mode == FAST_BINARY_SEARCH)
		{
		puts("Fast Binary Search");
		searcher = new fast_binary_search;
		}
	else if (mode == FAST_BINARY_SEARCH_AVX512)
		{
		exit(printf("Under development\n"));
//		searcher = new fast_binary_search_avx512;
		}
	else if (mode == HAMMING_DISTANCE)
		{
		puts("Hamming Distance\n");
		searcher = new hamming_distance;
		}

	std::cout << "Indexing\n";
	searcher->make_index(workload.genome_guides);
	std::cout << "Search\n";

	/*
		Start timing.
	*/
    auto start2 = std::chrono::steady_clock::now();

	/*
		Create space to put the results
	*/
	all_positions.resize(thread_count);

	/*
		Allocate the thread pool
	*/
	std::vector<std::thread> threads;

	/*
		Launch each thread
	*/
	std::cout << "Launching " << thread_count << " threads\n";
	for (size_t i = 0; i < thread_count; i++)
		threads.push_back(std::thread(&finder::process_chunk, searcher, std::ref(workload)));

	/*
		Wait for each thread to terminate
	*/
	for (auto &thread : threads)
		thread.join();

	/*
		Stop timing
	*/
	auto end2 = std::chrono::steady_clock::now(); // End timing
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count();

	delete searcher;

	/*
		Dump the statistics
	*/
	std::cout << "Number of test guides: " << TESTSIZE << '\n';
	std::cout << "Number of test guides with matches within 4: " << workload.hits << '\n';
	std::cout << "Best score: " << workload.best_score << " " << packer_2bit.unpack_20mer(workload.best_20mer) << '\n';
	std::cout << "Mean Execution time (getvar+intersect): " << duration2 / 1000000.001 << " seconds" << '\n';
	auto end_main = std::chrono::steady_clock::now(); // End timing
	auto duration_main = std::chrono::duration_cast<std::chrono::microseconds>(end_main - start_main).count();
	std::cout << "Total Execution time: " << duration_main / 1000000.001 << " seconds" << '\n';

	return 0;
	}
