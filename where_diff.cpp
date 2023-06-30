#include <bitset>
#include <iostream>
#include <algorithm>

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
	FOLDED_TO_LIST()
	----------------
*/
uint64_t folded_to_list(uint64_t *answer, uint64_t folded)
	{
	uint64_t bit_position[] = {0x000002, 0x000008, 0x000020, 0x000080, 0x000200, 0x000800, 0x002000, 0x008000, 0x020000, 0x080000, 0x000001, 0x000004, 0x000010, 0x000040, 0x000100, 0x000400, 0x001000, 0x004000, 0x010000, 0x040000};

	uint64_t *into = answer;
	for (uint64_t bit = 0; bit < 20; bit++)
		if (folded & bit_position[bit])
			*into++ = bit;

	return into - answer;
	}

/*
	UNFOLD()
	--------
	The bit positions are:
	9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
*/
uint64_t unfold(uint64_t folded)
	{
	uint64_t bit_position[] = {0x000002, 0x000008, 0x000020, 0x000080, 0x000200, 0x000800, 0x002000, 0x008000, 0x020000, 0x080000, 0x000001, 0x000004, 0x000010, 0x000040, 0x000100, 0x000400, 0x001000, 0x004000, 0x010000, 0x040000};
	uint64_t answer = 0;

	for (uint64_t bit = 0; bit < 20; bit++)
		if (folded & bit_position[bit])
			answer |= 1 << bit;

	return answer;
	}

/*
	FOLDED_WHERE_DIFFERENT()
	------------------------
*/
uint64_t folded_where_different(uint64_t guide, uint64_t variant)
	{
	/*
		Compute the "pairs" where there is a difference
	*/
	uint64_t xor_set = guide ^ variant;

	/*
		Take the high bits and the low bits and OR them resulting in only the high bit of the
		pair being set if either bit is set (i.e. there is a change).
	*/
	uint64_t diffs = (((xor_set & 0x5555555555555555) << 1) | (xor_set & 0xAAAAAAAAAAAAAAAA));

	/*
		Take the low 20 bits and fold into it the top 20 bits.  The resultant bitstring will be in the order:
		9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
	*/
	uint64_t folded = (diffs & 0xFFFFF) | (diffs >> 21);

	return folded;
	}

/*
	CALCULATE_MITLOCAL_SCORE()
	--------------------------
	Supplied by Dimitri Perrin

	This function calculates the local MIT score based on the positions of mismatches.
	https://dx.doi.org/10.1038/nbt.2647

	@param mismatch_array, An int array that contains the position of mismatches
	@param length, The length of the param `mismatch_array`
	@return The local MIT score
*/
double calculate_MITLocal_score(uint64_t *mismatch_array, uint64_t length)
	{
	uint64_t i;
	double T1 = 1.0, T2, T3, d = 0.0, score;
	/* Mismatch penalty array */
	double M[] = {0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

	/* 1st term */
	for (i = 0; i < length; ++i)
		T1 = T1 * (1.0 - M[mismatch_array[i]]);

	/* 2nd term */
	if (length == 1)
		d = 19.0;
	else
		{
		for (i = 0; i < length - 1; ++i)
			d += mismatch_array[i + 1] - mismatch_array[i];
		d = d / (length - 1);
		}
	T2 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1);

	/* 3rd term */
	T3 = 1.0 / (length * length);

	/* Local score */
	score = T1 * T2 * T3 * 100;

	return score;
	}

double MITLocal_score[1024 * 1024];
/*
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	kmer_encoding_table[(size_t)'A'] = 0;
	kmer_encoding_table[(size_t)'C'] = 1;
	kmer_encoding_table[(size_t)'G'] = 2;
	kmer_encoding_table[(size_t)'T'] = 3;

	uint64_t positions[20];
	double sum_of_scores = 0;
	for (uint64_t difflist = 1; difflist < 1024 * 1024; difflist++)
		{
		uint64_t length = folded_to_list(positions, difflist);
		double score = calculate_MITLocal_score(positions, length);
		MITLocal_score[difflist] = score;
		sum_of_scores += score;
		}

	std::cout << sum_of_scores << "\n";
//                               98765432109876543210
	uint64_t first =  pack20mer("TAAAAAAAAAAAAAAAAAAA");
	uint64_t second = pack20mer("AAAAAAAAAAAAAAAAAAAA");

	uint64_t diffs = folded_where_different(first, second);

	std::cout << std::bitset<20>(unfold(diffs)) << "\n";
	std::cout << std::bitset<20>(where_different(first, second)) << "\n";

	//
	// 39:37 of https://qut.zoom.us/rec/play/E_58mGQqRGqvxz0jCj9EtVib5pwc0i9hc1ov1V7vrbBeCxUwFn7iHOqcpef1vCWXBpC6liOUWXrPFYFA.f0hk00hXyISe0nHu?canPlayFromShare=true&from=share_recording_detail&continueMode=true&componentName=rec-play&originRequestUrl=https%3A%2F%2Fqut.zoom.us%2Frec%2Fshare%2F6_0N2sKRYUBwJ2rpmr51ol6eWlprGd-XVIXTb8J1_ff1PTz88mmfX-kqtDxz5B7p.JYjw2jy_lVvnMtRZ
	//
	// 35:03 of video
	// Global score = 100 / (100 + sum(MITLocal_scores))
	//
	std::cout << "Score:" << MITLocal_score[diffs] << "\n";

	return 0;
	}
