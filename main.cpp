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

#include "file.h"
#include "finder.h"
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
std::vector<uint64_t> read_guides(const std::string &filename, PACKER pack_20mer, bool map = true)
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
	auto new_end = std::unique(packed_guides.begin(), packed_guides.end());
	packed_guides.resize(std::distance(packed_guides.begin(), new_end));

	return packed_guides;
	}

/*
	USAGE()
	-------
*/
int usage(const char *exename)
	{
	std::cout << "Usage:" << exename << "[-b | -B | -h] [-f<filename>] [-t<threadcount>\n";
	std::cout << "       -b for binary search [default]\n";
	std::cout << "       -B for binary search using AVX512 instructions\n";
	std::cout << "       -h for hamming distance\n";
	std::cout << "       -t<threads> the number of threads to use when searching [default = corecount (inc hyperthreads)]\n";
	std::cout << "       -f<filename> use <filename> as the genome\n";

	return 0;
	}

/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
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
			else if (strncmp(argv[arg], "-t", 2) == 0)
				{
				if (strlen(argv[arg]) == 2)
					thread_count = atol(argv[++arg]);
				else
					thread_count = atol(argv[arg] + 2);
				}
			else
				{
				std::cout << "unknown paramter:" << argv[arg] << "\n";
				exit(usage(argv[0]));
				}
			}
		}
	else
		exit(usage(argv[0]));

	std::cout << "Processing using: " << mode_name[mode] << "\n";
	std::cout << "Using file      : " << guides_filename << "\n";

	/*
		Read the guides from the file.
		NOTE:read_guides() sorts the list after loading and before returning it.
	*/
//	std::string guides_filename = "OryzaSativaGuides.txt";
	auto time_io_start = std::chrono::steady_clock::now(); // Start timing
	std::vector<uint64_t> packed_genome_guides;
	if (mode == FAST_BINARY_SEARCH || mode == FAST_BINARY_SEARCH_AVX512)
		packed_genome_guides = read_guides(guides_filename, packer_2bit.pack_20mer);
	else
		packed_genome_guides = read_guides(guides_filename, packer_3bit.pack_20mer);
	auto time_io_end = std::chrono::steady_clock::now(); // Stop timing
	auto time_io_duration = std::chrono::duration_cast<std::chrono::microseconds>(time_io_end - time_io_start).count();
	std::cout << "Load encode time: " << time_io_duration / 1000000.001 << " seconds" << '\n';

	/*
		Generate some samples (sometimes random, and sometimes not, so keep both versions)
	*/
	std::vector<uint64_t> test_guides;
	if (true)
		{
		TESTSIZE = 10'000;
		test_guides.resize(TESTSIZE);
//		select_random_vectors(test_guides, packed_genome_guides);
		select_pseudo_random_vectors(test_guides, packed_genome_guides);
		sort(test_guides.begin(), test_guides.end());
		}
	else
		{
		TESTSIZE = packed_genome_guides.size();
		test_guides.resize(TESTSIZE);
		for (size_t which = 0; which < TESTSIZE; which++)
			test_guides[which] = packed_genome_guides[which];
		}

	std::cout << "Loaded " << test_guides.size() << " test guides, " <<  packed_genome_guides.size() << " genome guides" << '\n';

	/*
		Generate the index over the guides.
	*/
	finder *searcher;
	if (mode == FAST_BINARY_SEARCH)
		searcher = new fast_binary_search;
	else if (mode == FAST_BINARY_SEARCH_AVX512)
		{
		exit(printf("Under development\n"));
//		searcher = new fast_binary_search_avx512;
		}
	else // if (mode == HAMMING_DISTANCE)
		{
		exit(printf("Under development\n"));
//		searcher = new hamming_distance;
		}

	searcher->make_index(packed_genome_guides);

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
		Work out how big each thread-chunk should be
	*/
	const size_t chunk_size = test_guides.size() / thread_count;
	size_t start_index = 0;

	/*
		Launch each thread
	*/
	std::cout << "Launching " << thread_count << " threads\n";
	for (size_t i = 0; i < thread_count - 1; i++)
		{
		threads.push_back(std::thread(&finder::process_chunk, searcher, start_index, start_index + chunk_size, std::ref(test_guides), std::ref(packed_genome_guides), std::ref(all_positions[i])));
		start_index += chunk_size;
		}
	threads.push_back(std::thread(&finder::process_chunk, searcher, start_index, test_guides.size(), std::ref(test_guides), std::ref(packed_genome_guides), std::ref(all_positions[thread_count - 1])));

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
		Calculate the statistics
	*/
	uint64_t match_count = 0;
	double min_score = std::numeric_limits<double>::max();
	uint64_t best = 0;
	double max_score = std::numeric_limits<double>::min();
	for (int i = 0; i < all_positions.size(); i++)					// per thread
		for (int j = 0; j < all_positions[i].size(); j++)			// per guide (i.e. per guide per thread)
			{
			match_count++;
			double score = all_positions[i][j].score;
			if (score > max_score)
				max_score = score;
			if (score < min_score)
				{
				min_score = score;
				best = all_positions[i][j].sequence;
				}
			}

	/*
		Dump the statistics
	*/
	std::cout << "Number of test guides: " << TESTSIZE << '\n';
	std::cout << "Number of test guides with matches within 4: " << match_count << '\n';
	std::cout << "Min score: " << min_score << " " << packer_2bit.unpack_20mer(best) << '\n';
	std::cout << "Max score: " << max_score << '\n';
	std::cout << "Mean Execution time (getvar+intersect): " << duration2 / 1000000.001 << " seconds" << '\n';
	auto end_main = std::chrono::steady_clock::now(); // End timing
	auto duration_main = std::chrono::duration_cast<std::chrono::microseconds>(end_main - start_main).count();
	std::cout << "Total Execution time: " << duration_main / 1000000.001 << " seconds" << '\n';

	return 0;
	}
