#include <stdint.h>

#include <chrono>
#include <iostream>

	static uint64_t kmer_encoding_table[256];

	/*
		PACK20MER()
		-----------
	*/
	inline uint64_t pack20mer(const char *sequence)
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
	inline std::string unpack20mer(uint64_t packed_sequence)
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


uint64_t reverseComplement(uint64_t packedKmer, int k)
	{
	uint64_t reverseComplement = 0;
	uint64_t mask = 3; // Bit mask for extracting 2 bits at a time
	for (int i = 0; i < k; i++)
		{
		uint64_t base = packedKmer & mask;
		reverseComplement <<= 2;
		reverseComplement |= (~base & mask);
		packedKmer >>= 2;
		}
	// Reverse the order of the bases
	uint64_t reversedKmer = 0;
	while (reverseComplement != 0)
		{
		reversedKmer <<= 2;
		reversedKmer |= (reverseComplement & mask);
		reverseComplement >>= 2;
		}
	return reversedKmer;
	}

uint64_t at_reverse_complement(uint64_t kmer, size_t bases)
	{
	uint64_t result = 0;
	/*
		Compute the complement (A<->T, C<->G), which is a bit-fit because A=00, T=11, C=01, G=10
	*/
	uint64_t complement = ~kmer;

	/*
		Reverse the order of the bases
	*/
	for (uint64_t base = 0; base < bases; base++, complement >>= 2)
		result = (result << 2) | (complement & 0x03);

	return result;
	}

/*
	MAIN()
	------
*/
int main(void)
	{
	uint64_t checksum = 0;
	kmer_encoding_table[(size_t)'A'] = 0;
	kmer_encoding_table[(size_t)'C'] = 1;
	kmer_encoding_table[(size_t)'G'] = 2;
	kmer_encoding_table[(size_t)'T'] = 3;

	auto start = std::chrono::steady_clock::now(); // Start timing
	for (uint64_t times = 0; times < 100'000'000; times++)
		checksum += reverseComplement(times, 20);
	auto end = std::chrono::steady_clock::now(); // End timing

	std::cout << "checksum:" << checksum;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Total Execution time: " << duration / 1000000.001 << " seconds" << '\n';

	{
	uint64_t val = 0x1B;  // 0001 1011
	uint64_t rev_com = at_reverse_complement(val, 4);
	printf("%02llX %02llX\n", val, rev_com);
	}
	{
	uint64_t val = 0x2B13;  // 10 1011 0001 0011
	uint64_t rev_com = at_reverse_complement(val, 7);
	printf("%04llX %04llX\n", val, rev_com);
	}

	return 0;
	}
