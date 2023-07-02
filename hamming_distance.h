/*
	HAMMING_DISTANCE.H
	------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief The hamming distance method of finding variants of k-mers in a genome.
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/
#if defined(__APPLE__) || defined(__linux__)
	#define __popcnt64 __builtin_popcountll
#endif

/*
	CLASS HAMMING_DISTANCE
	----------------------
*/
/*!
	@brief The hamming distance method of finding variants of k-mers in a genome.
*/
class hamming_distance : public finder
	{
	private:
		/*
			HAMMING_DISTANCE::DISTANCE()
			----------------------------
		*/
		/*!
			@brief Compute the  hamming distance between two integers.
			@param a [in] The first integer.
			@param b [in] The second integer.
			@returns The number of bit positions in which they differ
		*/
		inline uint64_t distance(uint64_t a, uint64_t b)
			{
			return __popcnt64(a ^ b);
			}

		/*
			HAMMING_DISTANCE::COMPUTE_HAMMING_SET()
			---------------------------------------
		*/
		/*!
			@brief Run through the genome looking for matches.
			@param twice_the_max_distance [in] Twice the hamming distance looked for (will find all matches where the hamming distance <= twice\_the\_max\_distance / 2).
			@param key [in] The kmer to look for.
			@param data [in] The genome to look in.
			@param data\_length [in] The numner of k-mers in the genome being looked in.
			@param positions [out] The result set.
		*/
		void compute_hamming_set(uint64_t twice_the_max_distance, uint64_t key, const uint64_t *data, size_t data_length, std::vector<const uint64_t *> &positions)
			{
			const uint64_t *end = data + data_length;

			for (const uint64_t *which = data; which < end; which++)
				if (distance(key, *which) <= twice_the_max_distance)
					positions.push_back(which);
			}

	public:
		/*
			HAMMING_DISTANCE::~HAMMING_DISTANCE()
			-------------------------------------
		*/
		virtual ~hamming_distance()
			{
			/* Nothing */
			}

		/*
			HAMMING_DISTANCE::MAKE_INDEX()
			------------------------------
		*/
			/*!
				@brief Do any necessary indexing of the genome before starting the search process (for hamming distance, this does nothing) .
				@param genome [in] The genome being searched.
			*/

		virtual void make_index(const std::vector<uint64_t> &genome)
			{
			/* Nothing */
			}

		/*
			HAMMING_DISTANCE::PROCESS_CHUNK()
			---------------------------------
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
				compute_hamming_set(8, test_guides[which], packed_genome_guides.data(), packed_genome_guides.size(), answer[which]);
			}
	};
