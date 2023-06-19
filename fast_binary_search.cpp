/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#if defined(__APPLE__) || defined(__linux__)
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
#endif

#include <mutex>
#include <thread>
#include <chrono>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <algorithm>

namespace fast_binary_search
	{
	/*
		STRUCT START_END
		----------------
	*/
	typedef struct
		{
		const uint64_t *start;
		const uint64_t *end;
		} start_end;

	/*
		The number of bases to use as the key to the index.
	*/
	constexpr size_t index_width_in_bases = 14;
	constexpr size_t base_width_in_bits = 2;

	/*
		Fast lookup for decoding encoded kmers
	*/
	static uint64_t kmer_encoding_table[256];

	/*
		The index
	*/
	std::vector<start_end>index;

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
		index.resize((size_t)pow(4, index_width_in_bases));

		/*
			Get the end of each group by passing over the data once
		*/
		const uint64_t *end = data + length;
		for (const uint64_t *current = data; current < end; current++)
			index[*current >> ((20 - index_width_in_bases) * base_width_in_bits)].end = current;

		/*
			Get the start of each group by passing over the index and setting the start to one after the previosu end
		*/
		const uint64_t *previous = data;
		start_end *index_end = &index[index.size()];
		for (start_end *current = &index[0]; current < index_end; current++)
			{
			if (current->end != nullptr)
				{
				current->start = previous;
				previous = current->end + 1;
				}
			}
		}

	/*
		COMPUTE_INTERSECTION_LIST
		-------------------------
	*/
	void compute_intersection_list(const uint64_t *key, size_t key_length, const uint64_t *data, size_t data_length, std::vector<uint64_t> &matches, std::vector<size_t> &positions)
		{
		const uint64_t *key_end = key + key_length;
		for (const uint64_t *current_key = key; current_key < key_end; current_key++)
			{
			size_t index_key = *current_key >> ((20 - index_width_in_bases) * base_width_in_bits);
			if (index[index_key].start != nullptr)
				{

//std::cout << index[index_key].start << " - " << (index[index_key].end + 1);
//std::cout << " -> ";
//std::cout << index[index_key].start << " - " << index[index_key + 1].start;
//std::cout << "\n";

				const uint64_t *found = std::lower_bound(index[index_key].start, index[index_key].end + 1, *current_key);
				if (*found == *current_key)
					{
					matches.push_back(*found);
					positions.push_back(found - data);
					}
				}
			}
		}

	/*
		PACK20MER()
		-----------
	*/
	uint64_t pack20mer(const char *sequence)
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
	std::string unpack20mer(uint64_t packed_sequence)
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
					{
					contents = (char *)malloc(details.st_size);
					if (fread(contents, details.st_size, 1, fp) != 1)
						{
						free(contents);
						contents = NULL;
						}
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

		std::sort(packed_guides.begin(),packed_guides.end());
		free(data);

		return packed_guides;
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

			uint64_t was = (sequence >> shifter) & 3;
			uint64_t knock_out = (sequence & ~(3ULL << shifter));
			switch(was)
				{
				case 0:
					generate_variations_binary(knock_out | (1ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (2ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (3ULL << shifter), variations, replacements + 1, i + 1);
					break;
				case 1:
					generate_variations_binary(knock_out | (0ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (2ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (3ULL << shifter), variations, replacements + 1, i + 1);
					break;
				case 2:
					generate_variations_binary(knock_out | (0ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (1ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (3ULL << shifter), variations, replacements + 1, i + 1);
					break;
				case 3:
					generate_variations_binary(knock_out | (0ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (1ULL << shifter), variations, replacements + 1, i + 1);
					generate_variations_binary(knock_out | (2ULL << shifter), variations, replacements + 1, i + 1);
					break;
				}
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
	Mutexes
*/
std::mutex matches_mutex;
std::mutex positions_mutex;

/*
	Allocate space for the final set of results
*/
std::vector<std::vector<uint64_t>> all_matches;
std::vector<std::vector<size_t>> all_positions;

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

//		std::sort(variations.begin() + 1, variations.end());		// +1 because the first element is the test guide

		std::vector<uint64_t> matches;
		std::vector<size_t> positions;
		fast_binary_search::compute_intersection_list(variations.data() + 1, variations.size() - 1, packed_genome_guides.data(), packed_genome_guides.size(), matches, positions);

		// Lock the mutexes before modifying the shared vectors
		std::lock_guard<std::mutex> matches_lock(matches_mutex);
		std::lock_guard<std::mutex> positions_lock(positions_mutex);

		all_matches.push_back(matches);
		all_positions.push_back(positions);

//if (positions.size() == 0)
//	puts("Empty");
		}
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
		Generate some random samples.
	*/
	constexpr size_t TESTSIZE = 10000;
	std::vector<std::string> test_guides(TESTSIZE);
	fast_binary_search::select_random_vectors(test_guides, packed_genome_guides);

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
	for (int i = 0; i < thread_count - 1; i++)
		{
		threads.push_back(std::thread(process_chunk, start_index, start_index + chunk_size, std::ref(test_guides), std::ref(packed_genome_guides)));
//		process_chunk(start_index, start_index + chunk_size, test_guides, packed_genome_guides);
		start_index += chunk_size;
		}
	threads.push_back(std::thread(process_chunk, start_index, test_guides.size(), std::ref(test_guides), std::ref(packed_genome_guides)));
//	process_chunk(start_index, test_guides.size(), test_guides, packed_genome_guides);

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
	size_t total_count = all_matches.size();
	size_t empty_count = std::count_if(all_positions.begin(), all_positions.end(), [](const std::vector<size_t> &positions)
		{
		return positions.empty();
		});

	double sum_matches = 0.0;
	size_t min_matches = std::numeric_limits<size_t>::max();
	size_t max_matches = 0;

	for (const auto& matches : all_matches)
		{
		size_t count = matches.size();
		sum_matches += count;
		min_matches = std::min(min_matches, count);
		max_matches = std::max(max_matches, count);
		}
	double mean_matches = sum_matches / total_count;

	/*
		Calculate the standard deviation
	*/
	double variance = 0.0;
	for (const auto& matches : all_matches)
		{
		double diff = matches.size() - mean_matches;
		variance += diff * diff;
		}
	variance /= total_count;
	double std_dev_matches = std::sqrt(variance);

	/*
		Dump the statistics
	*/
	std::cout << "Mean number of matches: " << mean_matches << '\n';
	std::cout << "Standard deviation of matches: " << std_dev_matches << '\n';
	std::cout << "Minimum matches: " << min_matches << '\n';
	std::cout << "Maximum matches: " << max_matches << '\n';
	std::cout << "Number of unique test guides: " << empty_count << '\n';
	std::cout << "Number of test guides with off-targets: " << TESTSIZE - empty_count << '\n';
	std::cout << "Mean Execution time (getvar+intersect): " << duration2 / 1000000.001 << " seconds" << '\n';
	auto end_main = std::chrono::steady_clock::now(); // End timing
	auto duration_main = std::chrono::duration_cast<std::chrono::microseconds>(end_main - start_main).count();
	std::cout << "Total Execution time: " << duration_main / 1000000.001 << " seconds" << '\n';

	return 0;
	}

