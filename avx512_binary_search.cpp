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
#include <algorithm>

/*
	AVX_PRINT()
	-----------
*/
void avx_print(const char *str, uint64_t *array, __m512i value)
	{
	std::cout << str;
	uint64_t data[8];
	_mm512_storeu_epi64(&data[0], value);
	for (int x = 0; x < 2; x++)
		std::cout << '[' << data[x] << ' ' << array[data[x]] << ']';
	std::cout << '\n';
	}

void avx_print_ptr_and_data(const char *str, uint64_t *array, __m512i value)
	{
	std::cout << str;
	uint64_t *data[8];
	_mm512_storeu_epi64(&data[0], value);
	for (int x = 0; x < 2; x++)
		std::cout << "[" << data[x] << " " << data[x] - array << ' ' << *data[x] << "]";

	std::cout << '\n';
	}


/*
	AVX_BINARY_SEARCH()
	-------------------
	Uses indexes NOT pointers
*/
__m512i avx_binary_search(uint64_t *array, __m512i lower, __m512i upper, __m512i key)
	{
	__m512i one = _mm512_set1_epi64(1);
	__mmask8 not_finished = 0;

	while ((not_finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)
		{
		__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);

#ifdef NEVER
avx_print(".LOWER:", array, lower);
avx_print(".MIDDL:", array, middle);
avx_print(".UPPER:", array, upper);
std::cout << "\n";
#endif
		__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, array, sizeof(uint64_t)), key);
		lower = _mm512_mask_blend_epi64 (not_finished & results, lower, middle);
		upper = _mm512_mask_blend_epi64 (not_finished & ~results, upper, middle);
		}

	return upper;
	}


/*
	AVX_BINARY_SEARCH()
	-------------------
	Uses pointers rather than indexes
*/
__m512i avx_binary_search(__m512i lower, __m512i upper, __m512i key, uint64_t *data)
	{
	__m512i eight = _mm512_set1_epi64(8);
	__m512i high_bits = _mm512_set1_epi64(~(uint64_t)7);
	__mmask8 not_finished = 0;

	while ((not_finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, eight), upper)) != 0x00)			// eight bytes is 1 * sizeof(uint64_t)
		{
		__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);

#ifdef NEVER
avx_print_ptr_and_data("-LOWER:", data, lower);
avx_print_ptr_and_data("-MIDDL:", data, middle);
avx_print_ptr_and_data("-UPPER:", data, upper);
std::cout << "\n";
#endif
		middle = _mm512_and_epi64(middle, high_bits);



		__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, nullptr, 1), key);

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
	std::sort(data.begin(), data.end());

	for (size_t d = 0; d < 16'000'000; d++)
		search_for.push_back(distribution(random));

	std::cout << "64\n";
	uint64_t odds = 0;
	auto stopwatch_start = std::chrono::steady_clock::now();
	for (uint64_t x = 0; x < search_for.size(); x++)
		{
		auto found_at = std::lower_bound(&data[0], &data[data.size() - 1], search_for[x]);
//		if ((found_at - &data[0]) & 0x01)
//			odds++;
		if (*found_at & 0x01)
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
		auto key = _mm512_loadu_epi64 (&search_for[x]);
		auto found_at = avx_binary_search(&data[0], lower, upper, key);

		int64_t mem[8];
		_mm512_store_epi64 (&mem, found_at);
		for (uint64_t result_pos = 0; result_pos < 8; result_pos++)
			{
//			if (mem[result_pos] & 0x01)
//				odds++;
			if (data[mem[result_pos]] & 0x01)
				odds++;
//			std::cout << search_for[x + result_pos] << ':' << mem[result_pos] << '\n';
			}
//break;
		}
	stopwatch_end = std::chrono::steady_clock::now();
	stopwatch_duration = std::chrono::duration_cast<std::chrono::microseconds>(stopwatch_end - stopwatch_start).count();
	std::cout << "odds  : " << odds << '\n';
	std::cout << "512-bit:" << stopwatch_duration / 1000000.001 << " seconds" << '\n';



{
	std::cout << "512 pointers\n";
	odds = 0;
	stopwatch_start = std::chrono::steady_clock::now();
	__m512i lower = _mm512_set1_epi64((uint64_t)(&data[-1]));
	__m512i upper = _mm512_set1_epi64((uint64_t)(&data[0] + data.size() - 1));
	for (uint64_t x = 0; x < search_for.size(); x += 8)
		{
		auto key = _mm512_loadu_epi64 (&search_for[x]);
		auto found_at = avx_binary_search(lower, upper, key, &data[0]);

		int64_t *mem[8];
		_mm512_store_epi64 (mem, found_at);
		for (uint64_t result_pos = 0; result_pos < 8; result_pos++)
			{
			if (*mem[result_pos] & 0x01)
				odds++;
//			std::cout << search_for[x + result_pos] << ':' << mem[result_pos] << '\n';
			}
		}
	stopwatch_end = std::chrono::steady_clock::now();
	stopwatch_duration = std::chrono::duration_cast<std::chrono::microseconds>(stopwatch_end - stopwatch_start).count();
	std::cout << "odds  : " << odds << '\n';
	std::cout << "512-bit:" << stopwatch_duration / 1000000.001 << " seconds" << '\n';
}

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
