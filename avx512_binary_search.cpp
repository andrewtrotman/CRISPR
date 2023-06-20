
#include <immintrin.h>

#include <bitset>
#include <iostream>

/*
	BINARY_SEARCH()
	---------------
*/
int64_t binary_search(uint64_t *array, uint64_t elements, uint64_t key)
	{
	int64_t lower = -1;
	int64_t upper = elements - 1;

	while (lower + 1 != upper)
		{
		int64_t middle = (lower + upper) / 2;

	std::cout << "lower:" << lower << " middle:" << middle << " upper:" << upper << "\n";

		if (array[middle] < key)
			lower = middle;
		else
			upper = middle;
		}
	return upper;
	}

void print(__m512i lower, __m512i middle, __m512i upper)
	{
	int64_t lower_mem[8];
	int64_t middle_mem[8];
	int64_t upper_mem[8];

	_mm512_store_epi64 (&lower_mem, lower);
	_mm512_store_epi64 (&middle_mem, middle);
	_mm512_store_epi64 (&upper_mem, upper);

	std::cout << "lower:" << lower_mem[0] << " middle:" << middle_mem[0] << " upper:" << upper_mem[0] << "\n";
}

/*
	AVX_BINARY_SEARCH()
	-------------------
*/
int64_t avx_binary_search(uint64_t *array, int64_t elements, uint64_t the_key)
	{
	__mmask8 finished = 0;
	__m512i one = _mm512_set1_epi64(1);
	__m512i key = _mm512_set1_epi64(the_key);
	__m512i lower = _mm512_set1_epi64(-1);
	__m512i upper = _mm512_set1_epi64(elements - 1);

	while ((finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)					// while (lower + 1 != upper)
		{
		__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);			// middle = (lower + upper) / 2

		print(lower, middle, upper);

//		__mmask8 results = _mm512_cmplt_epi64_mask (_mm512_i64gather_epi64(middle, array, sizeof(uint64_t)), key); // array[middle] < key

		__mmask8 results = _mm512_mask_cmplt_epi64_mask(finished, _mm512_i64gather_epi64(middle, array, sizeof(uint64_t)), key);  // array[middle] < key

//{
//std::bitset<8> f_bits(finished);
//std::cout << "finished:" << f_bits << '\n';
//
//std::bitset<8> r_bits(results);
//std::cout << "results:" << r_bits << '\n';
//}

/*
		if (array[middle] < key)
			lower = middle;
		else
			upper = middle;
*/


		lower = _mm512_mask_blend_epi64 (results, lower, middle);		// if true update lower
		upper = _mm512_mask_blend_epi64 (results, middle, upper);		// else update upper
		}
	return upper;
	}



/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	std::cout << (1/2) << "\n";

	uint64_t array[] = {3, 4, 5, 6, 7, 8, 9};
	uint64_t search_for[] = {3, 4, 5, 6, 7, 8, 9, 10};
	int n = sizeof(array) / sizeof(array[0]);
	int64_t found_at = 0;

	for (uint64_t x = 0; x < sizeof(search_for); x++)
		{
		std::cout << "64\n";
		found_at = binary_search(array, n, search_for[x]);
		std::cout << x << ":" << found_at << "\n";
		}


		std::cout << "512\n";
		found_at = avx_binary_search(array, n, x);
		std::cout << x << ":" << found_at << "\n";
		std::cout << "\n";
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
