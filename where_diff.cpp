/*
	WHERE_DIFF.CPP
	--------------
*/
#include <bitset>
#include <iostream>
#include <algorithm>

#include "score_mit_local.h"

static uint64_t kmer_encoding_table[256];

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
	WHERE_DIFFERENT()
	-----------------
	Calculate the difference between the guide and the variant using XOR
	then there is a change-bit in either position 0 or 1 of the pair, so shift one to the right and OR.
	Except that high its will spill over so turn them off first by ANDing with 01 before shifting.
*/
uint64_t where_different(uint64_t guide, uint64_t variant)
	{
	uint64_t xor_set = guide ^ variant;
	uint64_t diffs = (((xor_set & 0x5555555555555555) << 1) | xor_set);

	uint32_t answer = 0;
	for (size_t pos = 0; pos < 32; pos++)
		answer |= (diffs & (3ULL << (pos * 2))) ? 1 << pos : 0;

	return answer;
	}


/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	score_mit_local scorer;
	kmer_encoding_table[(size_t)'A'] = 0;
	kmer_encoding_table[(size_t)'C'] = 1;
	kmer_encoding_table[(size_t)'G'] = 2;
	kmer_encoding_table[(size_t)'T'] = 3;

	uint64_t first =  pack20mer("TAAAAAAAAAAAAAAAAAAA");
	uint64_t second = pack20mer("AAAAAAAAAAAAAAAAAAAA");

	//
	// 39:37 of https://qut.zoom.us/rec/play/E_58mGQqRGqvxz0jCj9EtVib5pwc0i9hc1ov1V7vrbBeCxUwFn7iHOqcpef1vCWXBpC6liOUWXrPFYFA.f0hk00hXyISe0nHu?canPlayFromShare=true&from=share_recording_detail&continueMode=true&componentName=rec-play&originRequestUrl=https%3A%2F%2Fqut.zoom.us%2Frec%2Fshare%2F6_0N2sKRYUBwJ2rpmr51ol6eWlprGd-XVIXTb8J1_ff1PTz88mmfX-kqtDxz5B7p.JYjw2jy_lVvnMtRZ
	//
	// 35:03 of video
	// Global score = 100 / (100 + sum(MITLocal_scores))
	//
	std::cout << "Score:" << scorer.score(first, second) << "\n";

	return 0;
	}
