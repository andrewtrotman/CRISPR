
#include <immintrin.h>

#include <bitset>
#include <iostream>
#include <random>

/*
	AVX_BINARY_SEARCH()
	-------------------
*/
__m512i avx_binary_search(uint64_t *array, int64_t elements, __m512i key)
	{
	__mmask8 finished = 0;
	__m512i one = _mm512_set1_epi64(1);
	__m512i lower = _mm512_set1_epi64(-1);
	__m512i upper = _mm512_set1_epi64(elements - 1);

	while ((finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)
		{
		__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);
		__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, array, sizeof(uint64_t)), key);
		lower = _mm512_mask_blend_epi64 (finished & results, lower, middle);
		upper = _mm512_mask_blend_epi64 (finished & ~results, upper, middle);
		}
	return upper;
	}

/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	constexpr int seed = 17;
	std::mt19937 random(seed);
	std::uniform_int_distribution<int> distribution(0, guides.size() - 1);



	std::cout << (1/2) << "\n";

	uint64_t array[] = {3, 4, 5, 6, 7, 8, 9};
	uint64_t search_for[] = {0, 3, 4, 5, 6, 7, 8, 9};
	int n = sizeof(array) / sizeof(array[0]);
	

	for (uint64_t x = 0; x < sizeof(search_for) / sizeof(*search_for); x++)
		{
		std::cout << "64\n";
		auto found_at = binary_search(array, n, search_for[x]);
		std::cout << search_for[x] << "***FOUND:" << found_at << "\n";
		}

		{
		std::cout << "512\n";
		auto key = _mm512_load_epi64 (search_for);
		auto found_at = avx_binary_search(array, n, key);

		int64_t mem[8];

		_mm512_store_epi64 (&mem, found_at);
	
		for (uint64_t x = 0; x < sizeof(search_for) / sizeof(*search_for); x++)
			std::cout << search_for[x] << "***FOUND:" << mem[x] << "\n";
		}

	return 0;
	}

#ifdef NEVER
	/*
		The first occurance of key is at position upper if key is present in the array.
		Else, if upper > elements then we've gone off the top of the array.
	*/
	uint64_t answer = upper;
	if ((answer > elements) || (array[answer] != key))
		answer = -1;

	return answer;
#endif
