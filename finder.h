
/*
	FINDER.H
	--------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief Pure virtual base class for differnt finding techniques
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <vector>

#include "job.h"
#include "sequence_score_pair.h"

/*
	CLASS FINDER
	------------
*/
/*!
	@brief Pure virtual base class for differnt finding techniques
*/
class finder
	{
	public:
		/*
			FINDER::~FINDER()
			-----------------
		*/
		virtual ~finder()
			{
			/* Nothing */
			}
			
		/*
			FINDER::MAKE_INDEX()
			--------------------
		*/
			/*!
				@brief Do any necessary indexing of the genome before starting the search process (for hamming distance, this does nothing) .
				@param genome [in] The genome being searched.
			*/

		virtual void make_index(const std::vector<uint64_t> &genome) = 0;

		/*
			FINDER::PROCESS_CHUNK()
			-----------------------
		*/
		/*!
			@brief THREAD SAFE.  Search the genome for the given guides in test\_guides
			@param start [in] Where in the test_guides to start searching from.
			@param end [in] Where in the test_guides to start searching to.
			@param test_guides [in] The list of guides to look for - and we only look for those between start and end.
			@param packed_genome_guides [in] The genome to look in.
			@param answer [out] The result set
		*/
		virtual void process_chunk(job &workload, std::vector<uint64_t> &packed_genome_guides, std::vector<sequence_score_pair> &answer) = 0;
	};
