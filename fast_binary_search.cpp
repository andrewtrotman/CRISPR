/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include <math.h>
#include <stdio.h>
#include <stdint.h>

#if defined(__APPLE__) || defined(__linux__)
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
#endif

#include <vector>
#include <string>
#include <random>
#include <iostream>

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
	constexpr size_t index_width_in_bases = 5;
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
			index[*current >> ((19 - index_width_in_bases) * base_width_in_bits)].end = current;

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
				previous = current->end;
				}
			}
		}

	/*
		COMPUTE_INTERSECTION_LIST
		-------------------------
	*/
	void compute_intersection_list(const uint64_t *key, size_t key_length, const uint64_t *data, size_t data_length, std::vector<uint64_t> &matches, std::vector<size_t> &positions)
		{
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
		std::vector<uint64_t> packedGenomeGuides;
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
			packedGenomeGuides.push_back(pack20mer(guide));
			guide = strchr(guide, '\n');
			}
		while (guide != NULL && *(guide + 1) != '\0');

		std::sort(packedGenomeGuides.begin(),packedGenomeGuides.end());
		free(data);
		return packedGenomeGuides;
		}

	/*
		SELECT_RANDOM_VECTORS()
		-----------------------
	*/
	void select_random_vectors(std::vector<std::string> &selected, const std::vector<uint64_t>& guides)
	{
	constexpr int seed = 13;
	std::mt19937 random(seed);
	std::uniform_int_distribution<int> distribution(0, guides.size() - 1);

	for (size_t which = 0; which < selected.size(); which++)
		selected[which] = unpack20mer(guides[distribution(random)]);
	}
}

/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	fast_binary_search::kmer_encoding_table[(size_t)'A'] = 0;
	fast_binary_search::kmer_encoding_table[(size_t)'C'] = 1;
	fast_binary_search::kmer_encoding_table[(size_t)'G'] = 2;
	fast_binary_search::kmer_encoding_table[(size_t)'T'] = 3;

	/*
		Read the guides from the file
	*/
	std::string guides_filename = "OryzaSativaGuides.txt";
std::cout << "Load:" << guides_filename << '\n';
	std::vector<uint64_t> packed_genome_guides = fast_binary_search::load_guides(guides_filename);

	/*
		Generate some random samples
	*/
	constexpr size_t TESTSIZE = 10;
	std::vector<std::string> test_guides(TESTSIZE);
std::cout << "select random vectors" << '\n';
	fast_binary_search::select_random_vectors(test_guides, packed_genome_guides);

	for (const auto &guide : test_guides)
		{
		std::cout << guide << "->" << fast_binary_search::pack20mer(guide.c_str()) << '\n';
		if (guide != fast_binary_search::unpack20mer(fast_binary_search::pack20mer(guide.c_str())))
			std::cout << "ERROR\n";
		}

	return 0;
	}

