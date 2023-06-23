/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#if defined(__APPLE__) || defined(__linux__)
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>

	#define __popcnt16 __builtin_popcount
#endif

#include <thread>
#include <chrono>
#include <vector>
#include <string>
#include <random>
#include <bitset>
#include <iostream>
#include <algorithm>

#define __AVX512F__ 1
#undef __AVX512F__

#ifdef __APPLE__
        #define forceinline __attribute__((always_inline)) inline
#elif defined(__GNUC__)
        #define forceinline __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
        #define forceinline __forceinline
#else
        #define forceinline inline
#endif

namespace fast_binary_search
	{
	/*
		The number of bases to use as the key to the index.
	*/
	constexpr size_t index_width_in_bases = 14;
//	constexpr size_t index_width_in_bases = 1;
	constexpr size_t base_width_in_bits = 2;

	/*
		Fast lookup for decoding encoded kmers
	*/
	static uint64_t kmer_encoding_table[256];

	/*
		The index
	*/
	std::vector<const uint64_t *>index;

	/*
		Prototypes
	*/
	std::string unpack20mer(uint64_t packed_sequence);

	/*
		COMPUTE_INDEX()
		---------------
		data must be sorted
	*/
	void compute_index(const uint64_t *data, size_t length)
		{
		/*
			Allocate space for the index
		*/
		index.resize((size_t)pow(4, index_width_in_bases) + 1);		// +1 because we need a terminator on the end

		/*
			Get the end of each group by passing over the data once.  But first make sure the index starts at the beginning and finishes at the end.
		*/
		index[0] = data;
		index[index.size() - 1] = data + length - 1;
		const uint64_t *end = data + length;
		for (const uint64_t *current = data; current < end - 1; current++)
			index[(*current >> ((20 - index_width_in_bases) * base_width_in_bits)) + 1] = current + 1;

		/*
			Remove nullptrs by setting the start of the current range to the start of the previous range. Now,
		*/
		const uint64_t **index_end = &index[index.size() - 1];
		const uint64_t **previous = &index[0];
		for (const uint64_t **current = &index[0]; current < index_end; current++)
			{
			if (*current == nullptr)
				*current = *previous;
			previous = current;
			}
		}

	/*
		PACK20MER()
		-----------
	*/
	forceinline uint64_t pack20mer(const char *sequence)
		{
		uint64_t packed = 0;

		for (int pos = 0; pos < 20; pos++)
			packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];

		return packed;
		}

	/*
		UNPACK20MER()
		-------------
	*/
	forceinline std::string unpack20mer(uint64_t packed_sequence)
		{
		std::string sequence;
		for (int32_t pos = 19; pos >= 0; pos--)
			{
			switch ((packed_sequence >> (pos * 2)) & 3)
				{
				case 0:
					sequence += 'A';
					break;
				case 1:
					sequence += 'C';
					break;
				case 2:
					sequence += 'G';
					break;
				case 3:
					sequence += 'T';
					break;
				}
			}
		return sequence;
		}

#ifdef __AVX512F__
	/*
		AVX_BINARY_SEARCH()
		-------------------
	*/
	forceinline __m512i avx_binary_search(__m512i lower, __m512i upper, __m512i key)
		{
		__m512i one = _mm512_set1_epi64(1);
		lower = _mm512_srli_epi64(lower, 3);					// Convert from pointers to indexes
		lower = _mm512_sub_epi64(lower, one);					// Necessary because the binary search counts from 1 (not from 0)
		upper = _mm512_srli_epi64(upper, 3);					// Convert from pointers to indexes
		__mmask8 not_finished;

		/*
			Binary search
		*/
		while ((not_finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)
			{
			__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);
			__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, nullptr, sizeof(uint64_t)), key);
			lower = _mm512_mask_blend_epi64 (not_finished & results, lower, middle);
			upper = _mm512_mask_blend_epi64 (not_finished & ~results, upper, middle);
			}

		return _mm512_slli_epi64(upper, 3);		// convert back into pointers
		}
#endif

	/*
		COMPUTE_INTERSECTION_LIST()
		---------------------------
	*/
	void compute_intersection_list(const uint64_t *key, size_t key_length, const uint64_t *data, size_t data_length, std::vector<const uint64_t *> &positions)
		{
		positions.push_back(nullptr);				// push the guide before working on the variants

		const uint64_t *key_end = key + key_length;
		const uint64_t *current_key = key + 1;

		/*
			Do the AVX512 stuff first
		*/
#ifdef __AVX512F__
		__m512i one = _mm512_set1_epi64(1);

		while (current_key + 8 < key_end)
			{
			/*
				Get the current key set (that we're looking for)
			*/
			__m512i key_set = _mm512_loadu_epi64(current_key);
			__m512i index_key_set = _mm512_srli_epi64(key_set, (20 - index_width_in_bases) * base_width_in_bits);

			/*
				Compute the start and end of the ranges
			*/
			__m512i start_set = _mm512_i64gather_epi64(index_key_set, &index[0], sizeof(uint64_t *));
			__m512i end_set = _mm512_i64gather_epi64(_mm512_add_epi64(index_key_set, one), &index[0], sizeof(uint64_t *));

			/*
				Do the AVX512-parallel binary search given the pointers
			*/
			__m512i found_pointers = avx_binary_search(start_set, end_set, key_set);

			/*
				At this point found_set is a set of pointers to the data - which may or may not match the key (same as the result of std::lower_bound())
			*/
			__m512i found_values = _mm512_i64gather_epi64(found_pointers, nullptr, 1);		// 1 because found_set is a set of pointers to 64-bit integers
			__mmask8 found_masks = _mm512_cmpeq_epi64_mask(key_set, found_values);

			if (found_masks != 0)
				{
				/*
					I really want to write directly into the positions vector, but it isn't clear how large it needs to be pre-allocated
				*/
				alignas(64) uint64_t *position[8];
				_mm512_mask_compressstoreu_epi64 (&position[0], found_masks, found_pointers);

				switch(__popcnt16(found_masks))
					{
					case 8:
						positions.push_back(position[7]);			// FALL THROUGH
					case 7:
						positions.push_back(position[6]);			// FALL THROUGH
					case 6:
						positions.push_back(position[5]);			// FALL THROUGH
					case 5:
						positions.push_back(position[4]);			// FALL THROUGH
					case 4:
						positions.push_back(position[3]);			// FALL THROUGH
					case 3:
						positions.push_back(position[2]);			// FALL THROUGH
					case 2:
						positions.push_back(position[1]);			// FALL THROUGH
					case 1:
						positions.push_back(position[0]);			// FALL THROUGH
					}
				}
			current_key += 8;
			}
#endif

		/*
			Then fall through to process the last few one at a time (if no AVX512 then do this for all guides)
		*/
		while (current_key < key_end)
			{
			size_t index_key = *current_key >> ((20 - index_width_in_bases) * base_width_in_bits);
			size_t mers_to_search = index[index_key + 1] - index[index_key];
			if (mers_to_search == 0)
				{
				/*
					Nothing
				*/
				}
//#ifdef NEVER
			else if (mers_to_search < 4)
				{
				/*
					If the length of the list is short then
					Linear search in an ordered list
				*/
				const uint64_t *where = index[index_key];
				while (*where < *current_key)
					where++;
				if (*where == *current_key)
					positions.push_back(where);
				}
//#endif
			else
				{
				/*
					If the length of the list is long then
					Binary search in an ordered list
				*/
				const uint64_t *found = std::lower_bound(index[index_key], index[index_key + 1], *current_key);
				if (*found == *current_key)
					positions.push_back(found);
				}
			current_key++;
			}
		}

	/*
		READ_ENTIRE_FILE()
		------------------
	*/
	char *read_entire_file(const char *filename)
		{
		FILE *fp;
		struct stat details;
		char *contents = NULL;

		if ((fp = fopen(filename, "rb")) != NULL)
			{
			if (fstat(fileno(fp), &details) == 0)
				if (details.st_size != 0)
					if ((contents = (char *)malloc(details.st_size + 1)) != NULL)
						{
						if (fread(contents, details.st_size, 1, fp) != 1)
							{
							free(contents);
							contents = NULL;
							}
						contents[details.st_size] = '\0';
						}
			fclose(fp);
			}
	return contents;
	}

	/*
		LOAD_GUIDES()
		-------------
	*/
	std::vector<uint64_t> load_guides(const std::string &filename)
		{
		std::vector<uint64_t> packed_guides;
		char *guide;
		char *data = read_entire_file(filename.c_str());
		if (data == NULL)
			{
			std::cerr << "Error opening guide file: " << filename << std::endl;
			exit(1);
			}

		guide = data - 1;
		do
			{
			guide++;
			packed_guides.push_back(pack20mer(guide));
			guide = strchr(guide, '\n');
			}
		while (guide != NULL && *(guide + 1) != '\0');

		std::sort(packed_guides.begin(), packed_guides.end());
		free(data);

		return packed_guides;
		}

	/*
		SELECT_PSEUDO_RANDOM_VECTORS()
		------------------------------
	*/
	void select_pseudo_random_vectors(std::vector<std::string> &selected, const std::vector<uint64_t> &guides)
		{
		uint64_t seed = 17;
		uint64_t m = guides.size() - 1;		// Modulus parameter
		uint64_t a = 3;						// Multiplier term
		uint64_t c = 3;						// Increment term
		uint64_t next = seed;

		for (uint64_t which = 0; which < selected.size(); which++)
			{
			selected[which] = unpack20mer(guides[next]);
			next = ((next * a) + c) % m;
			}
		}

	/*
		SELECT_RANDOM_VECTORS()
		-----------------------
	*/
	void select_random_vectors(std::vector<std::string> &selected, const std::vector<uint64_t>& guides)
		{
//		constexpr int seed = 13;
		constexpr int seed = 17;
		std::mt19937 random(seed);
		std::uniform_int_distribution<int> distribution(0, guides.size() - 1);

		for (size_t which = 0; which < selected.size(); which++)
			selected[which] = unpack20mer(guides[distribution(random)]);
		}

	/*
		GENERATE_VARIATIONS_BINARY()
		----------------------------
	*/
	void generate_variations_binary(uint64_t sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0)
		{
		variations.push_back(sequence);
		if (replacements == 4)
			return;
		for (uint64_t i = position; i < 20; i++)
			{
			uint64_t shifter = (19 - i) * 2;

			/*
				This optimisation is from Timothy Chappell who points out that by XORing the
				current base with 01, 02, and 03 you automatically get all the variants except the
				original sequence (but not necessarily in a predictable order).  So we don't need to
				mask and replace, we simply xor with the original sequence and call outselves recursively.
			*/
			generate_variations_binary(sequence ^ (1ULL << shifter), variations, replacements + 1, i + 1);
			generate_variations_binary(sequence ^ (2ULL << shifter), variations, replacements + 1, i + 1);
			generate_variations_binary(sequence ^ (3ULL << shifter), variations, replacements + 1, i + 1);
			}
		}

	/*
		GENERATE_VARIATIONS()
		---------------------
	*/
	void generate_variations(std::string& sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0)
		{
		generate_variations_binary(pack20mer(sequence.c_str()), variations, replacements, position);
		}
}

/*
	Allocate space for the final set of results
*/
std::vector<std::vector<const uint64_t *>> all_positions;

/*
	Command line parameters
*/
size_t TESTSIZE = 1000;

/*
	PROCESS_CHUNK()
	---------------
*/
void process_chunk(size_t start, size_t end, std::vector<std::string> &test_guides, std::vector<uint64_t> &packed_genome_guides)
	{
	for (size_t i = start; i < end; ++i)
		{
		std::vector<uint64_t> variations;
		fast_binary_search::generate_variations(test_guides[i], variations);
		fast_binary_search::compute_intersection_list(variations.data(), variations.size(), packed_genome_guides.data(), packed_genome_guides.size(), all_positions[i]);
		}
	}

/*
	USAGE()
	-------
*/
void usage(const char *exename)
	{
	std::cout << "Usage:" << exename << "[-b<index_width_in_bases>] [-t<TESTSIZE>]\n";
	}

/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	auto start_main = std::chrono::steady_clock::now(); // Start timing

	fast_binary_search::kmer_encoding_table[(size_t)'A'] = 0;
	fast_binary_search::kmer_encoding_table[(size_t)'C'] = 1;
	fast_binary_search::kmer_encoding_table[(size_t)'G'] = 2;
	fast_binary_search::kmer_encoding_table[(size_t)'T'] = 3;

	/*
		Read the guides from the file.
		NOTE:load_guides() sorts the list after loading and before returning it.
	*/
	std::string guides_filename = "OryzaSativaGuides.txt";
	std::vector<uint64_t> packed_genome_guides = fast_binary_search::load_guides(guides_filename);

	/*
		Generate some samples (sometimes random, and sometimes not, so keep both versions)
	*/
	std::vector<std::string> test_guides;
	if (true)
		{
		TESTSIZE = 1000;
		test_guides.resize(TESTSIZE);
		fast_binary_search::select_random_vectors(test_guides, packed_genome_guides);
		fast_binary_search::select_pseudo_random_vectors(test_guides, packed_genome_guides);
		}
	else
		{
		TESTSIZE = packed_genome_guides.size();
		test_guides.resize(TESTSIZE);
		for (size_t which = 0; which < TESTSIZE; which++)
			test_guides[which] = fast_binary_search::unpack20mer(packed_genome_guides[which]);
		}

	std::cout << "Loaded " << test_guides.size() << " test guides, " <<  packed_genome_guides.size() << " genome guides" << '\n';

	/*
		Generate the index over the guides.
	*/
	fast_binary_search::compute_index(&packed_genome_guides[0], packed_genome_guides.size());

	/*
		Start timing.
	*/
    auto start2 = std::chrono::steady_clock::now();

	/*
		Create space to put the results
	*/
	all_positions.resize(test_guides.size());

	/*
		Allocate the thread pool
	*/
	size_t thread_count = std::thread::hardware_concurrency();
	thread_count = 1;
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
		threads.push_back(std::thread(process_chunk, start_index, start_index + chunk_size, std::ref(test_guides), std::ref(packed_genome_guides)));
		start_index += chunk_size;
		}
	threads.push_back(std::thread(process_chunk, start_index, test_guides.size(), std::ref(test_guides), std::ref(packed_genome_guides)));

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

	/*
		Calculate the statistics
	*/
	double sum_matches = 0.0;
	size_t min_matches = all_positions[0].size() + 1;
	size_t max_matches = all_positions[0].size() - 1;
//	size_t total_count = 0;
	size_t empty_count = 0;
	for (int i = 0; i < all_positions.size(); i++)
		{
		size_t num_matches = all_positions[i].size() - 1;
		sum_matches += num_matches;
		if (num_matches > 0)
			empty_count++;
		if (num_matches > max_matches)
			max_matches = num_matches;
		if (num_matches < min_matches)
			min_matches = num_matches;
		if ((min_matches == max_matches) && i > 0)
			std::cout << "OOPS! min = " << min_matches << "max = " << max_matches << ", i = " << i << std::endl;
		}
//	total_count = all_positions.size() - empty_count;
	double mean_matches = sum_matches / all_positions.size();

	/*
		Dump the statistics
	*/
	std::cout << "Mean number of matches: " << mean_matches << '\n';
	std::cout << "Minimum matches: " << min_matches << '\n';
	std::cout << "Maximum matches: " << max_matches << '\n';
	std::cout << "Number of unique test guides: " << TESTSIZE - empty_count << '\n';
	std::cout << "Number of test guides with off-targets: " << empty_count << '\n';
	std::cout << "Mean Execution time (getvar+intersect): " << duration2 / 1000000.001 << " seconds" << '\n';
	auto end_main = std::chrono::steady_clock::now(); // End timing
	auto duration_main = std::chrono::duration_cast<std::chrono::microseconds>(end_main - start_main).count();
	std::cout << "Total Execution time: " << duration_main / 1000000.001 << " seconds" << '\n';

	return 0;
	}
