/*
	FAST_BINARY_SEARCH.H
	--------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief The fast indexed binary search method of finding variants of k-mers in a genome.
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <stdlib.h>
#include <stdint.h>

#include <vector>

#include "finder.h"

/*
	CLASS FAST_BINARY_SEARCH
	------------------------
*/
/*!
	@brief The fast indexed binary search method of matching.
*/
class fast_binary_search : public finder
	{
	protected:
		static constexpr size_t index_width_in_bases = 14;			// bases to use as the key in the index
		static constexpr uint64_t in_genome_reverse_bitmap_width_in_bases = 14;
		static constexpr uint64_t in_genome_reverse_bitmap_mask = (0x1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) - 1;
		static constexpr size_t base_width_in_bits = 2;			// must be 2 (2 bits per base in the uint64_t encoding)

		/*
			Fast lookup for encoding encoded kmers, and the index too.
		*/
		std::vector<uint64_t> in_genome_bitmap; 		 // Bitmap array with 2^32 bits (i.e. 16 bases wide)
		std::vector<uint64_t> in_genome_reverse_bitmap;  // Bitmap array whose size id dependant on in_genome_reverse_bitmap_width_in_bases
		std::vector<const uint64_t *>index;

	protected:
		/*
			FAST_BINARY_SEARCH::IN_GENOME_HEAD()
			------------------------------------
		*/
		/*!
			@brief Is the prefix of the k-mer anywhere in the genome being searched?
			@param key [in] The k-mer being checked.
			@returns true if the prefix is present, false otherwise.
		*/
		inline bool in_genome_head(uint64_t key) const
			{
			return in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64));
			}

		/*
			FAST_BINARY_SEARCH::IN_GENOME()
			-------------------------------
		*/
		/*!
			@brief Are both  the prefix of the k-mer and the suffix of the k-mer present in the genome?
			@param key [in] The k-mer being checked.
			@returns true if the prefix and suffix are both present, false otherwise.
		*/
		inline bool in_genome(uint64_t key) const
			{
			return
				(in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64)))
				&&
				(in_genome_reverse_bitmap[(key & in_genome_reverse_bitmap_mask) / 64] & (1ULL << ((key & in_genome_reverse_bitmap_mask) % 64)))
				;
			}

		/*
			FAST_BINARY_SEARCH::COMPUTE_INTERSECTION()
			------------------------------------------
		*/
		/*!
			@brief Compute the intersecton between a key and the genome the has already been indexed.
			@param key [in] the 20-mers being searched for.
			@param positions [out] The vector of match locations (pointers).
		*/
		void compute_intersection(const uint64_t key, std::vector<const uint64_t *> &positions) const;

		/*
			FAST_BINARY_SEARCH::COMPUTE_INTERSECTION_LIST()
			-----------------------------------------------
		*/
		/*!
			@brief Compute the intersecton between the keys and the genome the has already been indexed.
			@param key [in] the list of 20-mers being searched for.
			@param key_length [in] the number of k-mers encoded in key.
			@param positions [out] The vector of match locations (pointers).
		*/
		void compute_intersection_list(const uint64_t *key, size_t key_length, std::vector<const uint64_t *> &positions) const;

		/*
			FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
			-----------------------------------------
		*/
		/*!
			@brief Generate all 20-mers that are 0 to 4 base mutations away from sequence (and also "likely" to be in the genome begine searched)
			@param sequence [in] The sequence being mutated.
			@param variations [out] The vector of variations.
			@param replacements [in] The number of mutations that have already been applied to the sequence.
			@param position [in] The position to start mautating from
		*/
		void generate_variations(uint64_t sequence, std::vector<uint64_t> &variations, int replacements, int position) const;

		/*
			FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
			-----------------------------------------
		*/
		/*!
			@brief Generate all 20-mers that are 0 to 4 base mutations away from sequence (and also "likely" to be in the genome begine searched)
			@param sequence [in] The sequence being mutated.
			@param variations [out] The vector of variations.
		*/
		void generate_variations(uint64_t sequence, std::vector<uint64_t> &variations) const
			{
			generate_variations(sequence, variations, 0, 0);
			}

	public:
		/*
			FAST_BINARY_SEARCH::FAST_BINARY_SEARCH()
			----------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		fast_binary_search() :
			in_genome_bitmap((1ULL << 32) / 64, 0),
			in_genome_reverse_bitmap((1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) / 64, 0)
			{
			/*
				Nothing
			*/
			}

		/*
			FAST_BINARY_SEARCH::~FAST_BINARY_SEARCH()
			-----------------------------------------
		*/
		/*!
			@brief Destructor
		*/
		virtual ~fast_binary_search()
			{
			/* Nothing */
			}


		/*
			FAST_BINARY_SEARCH::MAKE_INDEX()
			--------------------------------
			data must be sorted
		*/
		/*!
			@brief Do any necessary indexing of the genome before starting the search process.
			@param genome [in] The genome being searched.
		*/
		virtual void make_index(const std::vector<uint64_t> &genome);

		/*
			FAST_BINARY_SEARCH::PROCESS_CHUNK()
			-----------------------------------
		*/
		/*!
			@brief THREAD SAFE.  Search the genome for the given guides in test\_guides
			@param start [in] Where in the test_guides to start searching from.
			@param end [in] Where in the test_guides to start searching to.
			@param test_guides [in] The list of guides to look for - and we only look for those between start and end.
			@param packed_genome_guides [in] The genome to look in.
			@param answer [out] The result set
		*/
		virtual void process_chunk(size_t start, size_t end, std::vector<uint64_t> &test_guides, std::vector<uint64_t> &packed_genome_guides, std::vector<std::vector<const uint64_t *>> &answer)
			{
			std::vector<uint64_t> variations;

			for (size_t which = start; which < end; which++)
				{
				variations.clear();
				generate_variations(test_guides[which], variations);
				compute_intersection_list(variations.data(), variations.size(), answer[which]);
				}
			}
	};
