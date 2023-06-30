#include <bitset>
#include <iostream>

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
	UNFOLD()
	--------
	9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
*/
uint64_t unfold(uint64_t folded)
	{
	uint64_t answer = 0;

	answer |= (((folded & 0x000002) >= 1 ? 1 : 0) << 0);
	answer |= (((folded & 0x000008) >= 1 ? 1 : 0) << 1);
	answer |= (((folded & 0x000020) >= 1 ? 1 : 0) << 2);
	answer |= (((folded & 0x000080) >= 1 ? 1 : 0) << 3);
	answer |= (((folded & 0x000200) >= 1 ? 1 : 0) << 4);
	answer |= (((folded & 0x000800) >= 1 ? 1 : 0) << 5);
	answer |= (((folded & 0x002000) >= 1 ? 1 : 0) << 6);
	answer |= (((folded & 0x008000) >= 1 ? 1 : 0) << 7);
	answer |= (((folded & 0x020000) >= 1 ? 1 : 0) << 8);
	answer |= (((folded & 0x080000) >= 1 ? 1 : 0) << 9);
	answer |= (((folded & 0x000001) >= 1 ? 1 : 0) << 10);
	answer |= (((folded & 0x000004) >= 1 ? 1 : 0) << 11);
	answer |= (((folded & 0x000010) >= 1 ? 1 : 0) << 12);
	answer |= (((folded & 0x000040) >= 1 ? 1 : 0) << 13);
	answer |= (((folded & 0x000100) >= 1 ? 1 : 0) << 14);
	answer |= (((folded & 0x000400) >= 1 ? 1 : 0) << 15);
	answer |= (((folded & 0x001000) >= 1 ? 1 : 0) << 16);
	answer |= (((folded & 0x004000) >= 1 ? 1 : 0) << 17);
	answer |= (((folded & 0x010000) >= 1 ? 1 : 0) << 18);
	answer |= (((folded & 0x040000) >= 1 ? 1 : 0) << 19);

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
	MAIN()
	------
*/
int main(int argc, const char *argv[])
	{
	kmer_encoding_table[(size_t)'A'] = 0;
	kmer_encoding_table[(size_t)'C'] = 1;
	kmer_encoding_table[(size_t)'G'] = 2;
	kmer_encoding_table[(size_t)'T'] = 3;

//                               98765432109876543210
	uint64_t first =  pack20mer("TTTTTTTTTTAAAAAAAAAT");
	uint64_t second = pack20mer("AAAAAAAAAAAAAAAAAATA");

	uint64_t diffs = folded_where_different(first, second);

	std::cout << std::bitset<20>(unfold(diffs)) << "\n";
	std::cout << std::bitset<20>(where_different(first, second)) << "\n";
	}
