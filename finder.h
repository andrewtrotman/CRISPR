
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
			@param workload [in] the queries, data, and frequencies
			@param answer [out] The result set
		*/
		virtual void process_chunk(job &workload) = 0;
	};
