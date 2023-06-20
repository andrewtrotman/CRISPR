/*
	AVX512_BINARY_SEARCH()
	----------------------
	Copyright (c) 2023 Andrew Trotman.
	Perform 8 binary seaches in parallel, one in each of the 8 64-bit parts of an AVX512 register.
*/
#include <immintrin.h>

#include <bitset>
#include <chrono>
#include <random>
#include <iostream>

/*
	AVX_BINARY_SEARCH()
	-------------------
*/
__m512i avx_binary_search(uint64_t *array, __m512i lower, __m512i upper, __m512i key)
	{
	__m512i one = _mm512_set1_epi64(1);
	__mmask8 not_finished = 0;

	while ((not_finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)
		{
		__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);
		__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, array, sizeof(uint64_t)), key);
		lower = _mm512_mask_blend_epi64 (not_finished & results, lower, middle);
		upper = _mm512_mask_blend_epi64 (not_finished & ~results, upper, middle);
		}

	return upper;
	}

/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	constexpr int data_size = 1'000'000;
	constexpr int seed = 17;
	std::mt19937 random(seed);
	std::uniform_int_distribution<int> distribution(0, data_size - 1);

	std::vector<uint64_t> data;
	std::vector<uint64_t> search_for;

	for (size_t d = 0; d < 100; d++)
		data.push_back(distribution(random));
	sort(data.begin(), data.end());

	for (size_t d = 0; d < 16'000'000; d++)
		search_for.push_back(distribution(random));

	std::cout << "64\n";
	uint64_t odds = 0;
	auto stopwatch_start = std::chrono::steady_clock::now();
	for (uint64_t x = 0; x < search_for.size(); x++)
		{
		auto found_at = std::lower_bound(&data[0], &data[data.size() - 1], search_for[x]);
		if ((found_at - &data[0]) & 0x01)
			odds++;
//		std::cout << search_for[x] << ':' << found_at - data.begin() << '\n';
		}
	auto stopwatch_end = std::chrono::steady_clock::now();
	auto stopwatch_duration = std::chrono::duration_cast<std::chrono::microseconds>(stopwatch_end - stopwatch_start).count();
	std::cout << "odds  : " << odds << '\n';
	std::cout << "64-bit: " << stopwatch_duration / 1000000.001 << " seconds" << '\n';

	std::cout << "512\n";
	odds = 0;
	stopwatch_start = std::chrono::steady_clock::now();
	__m512i lower = _mm512_set1_epi64(-1);
	__m512i upper = _mm512_set1_epi64(data.size() - 1);
	for (uint64_t x = 0; x < search_for.size(); x += 8)
		{
		auto key = _mm512_load_epi64 (&search_for[x]);
		auto found_at = avx_binary_search(&data[0], lower, upper, key);

		int64_t mem[8];
		_mm512_store_epi64 (&mem, found_at);
		for (uint64_t result_pos = 0; result_pos < 8; result_pos++)
			{
			if (mem[result_pos] & 0x01)
				odds++;
//			std::cout << search_for[x + result_pos] << ':' << mem[result_pos] << '\n';
			}
		}
	stopwatch_end = std::chrono::steady_clock::now();
	stopwatch_duration = std::chrono::duration_cast<std::chrono::microseconds>(stopwatch_end - stopwatch_start).count();
	std::cout << "odds  : " << odds << '\n';
	std::cout << "512-bit:" << stopwatch_duration / 1000000.001 << " seconds" << '\n';

	return 0;
	}

#ifdef NEVER
				const uint64_t *found = std::lower_bound(index[index_key], index[index_key + 1], *current_key);

	/*
		The first occurance of key is at position upper if key is present in the array.
		Else, if upper > elements then we've gone off the top of the array.
	*/
	uint64_t answer = upper;
	if ((answer > elements) || (array[answer] != key))
		answer = -1;

	return answer;
#endif
