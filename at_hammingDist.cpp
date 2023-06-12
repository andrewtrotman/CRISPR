/*
	FIX: 	line_count -= line_count % 4;		// so that we have a full number of SIMD registers
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cstdint>

#ifndef _MSC_VER
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>

	#define __popcnt64 __builtin_popcountll
#else
	#define __m256i_u __m256i
	#define __AVX512F__ 1
#endif

#include <immintrin.h>

#include "timer.h"


#define forceinline
/*
#ifdef __APPLE__
        #define forceinline __attribute__((always_inline)) inline
#elif defined(__GNUC__)
        #define forceinline __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
        #define forceinline __forceinline
#else
        #define forceinline inline
#endif
*/

char *read_entire_file(const char *filename)
	{
	FILE *fp;
	struct stat details;
	char *contents = NULL;

	if ((fp = fopen(filename, "rb")) != NULL)
		{
		if (fstat(fileno(fp), &details) == 0)
			if (details.st_size != 0)
				{
				contents = (char *)malloc(details.st_size);
				if (fread(contents, details.st_size, 1, fp) != 1)
					{
					free(contents);
					contents = NULL;
					}
				}
		fclose(fp);
		}
return contents;
}

char **buffer_to_list(char *buffer, uint64_t *lines)
{
char *pos, **line_list, **current_line;
uint64_t n_frequency, r_frequency;

n_frequency = r_frequency = 0;
for (pos = buffer; *pos != '\0'; pos++)
	if (*pos == '\n')
		n_frequency++;
	else if (*pos == '\r')
		r_frequency++;

*lines = r_frequency > n_frequency ? r_frequency : n_frequency;
current_line = line_list = new (std::nothrow) char * [(size_t)(*lines + 2)]; 		// +1 in case the last line has no \n; +1 for a NULL at the end of the list

if (line_list == NULL)		// out of memory!
	return NULL;

*current_line++ = pos = buffer;
while (*pos != '\0')
	{
	if (*pos == '\n' || *pos == '\r')
		{
		*pos++ = '\0';
		while (*pos == '\n' || *pos == '\r')
			pos++;
		*current_line++ = pos;
		}
	else
		pos++;
	}
/*
	We have a nasty case here.  If the last line has no CR/LF then we need to include it
	but shove a NULL on the next line, but if the last line has a CR/LF then we need to avoid
	adding a blank line to the end of the list.
	NOTE: its 2012 and its a bit late to be finding this bug!!!
*/
if (**(current_line - 1) == '\0')
	*(current_line - 1) = NULL;
*current_line = NULL;

*lines = current_line - line_list - 1;		// the true number of lines

return line_list;
}



forceinline uint64_t hammingDistance(uint64_t a, uint64_t b)
	{
	uint64_t xorResult = a ^ b;
//	return std::bitset<64>(xorResult).count();
	return __popcnt64(xorResult);
	}

void filterAndSortByHammingDistance(uint64_t a, const std::vector<uint64_t>& b, uint64_t maxDistance, std::vector<uint64_t>& filteredB) {
    // calculate the distance between a and every element in b, and return elements where dist(a,b)<maxDistance 
    filteredB.reserve(b.size());

    for (size_t i = 0; i < b.size(); i++) {
        uint64_t distance = hammingDistance(a, b[i]);
        if (distance <= maxDistance) {
            filteredB.push_back(b[i]);
        }
    }

    std::sort(filteredB.begin(), filteredB.end());
}

#ifdef __AVX512F__
/*
forceinline __m256i hs_popcount(const __m256i v)
	{
	const __m256i m1 = _mm256_set1_epi8(0x55);
	const __m256i m2 = _mm256_set1_epi8(0x33);
	const __m256i m4 = _mm256_set1_epi8(0x0F);

	const __m256i t1 = _mm256_sub_epi8(v,       (_mm256_srli_epi16(v,  1) &m1));
	const __m256i t2 = _mm256_add_epi8(t1 & m2, (_mm256_srli_epi16(t1, 2) &m2));
	const __m256i t3 = _mm256_add_epi8(t2, _mm256_srli_epi16(t2, 4)) &m4;
	return _mm256_sad_epu8(t3, _mm256_setzero_si256());
	}
*/

const __m256i lookup = _mm256_setr_epi8(
            /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
            /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
            /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
            /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,
            /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
            /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
            /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
            /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4
    );
const __m256i low_mask = _mm256_set1_epi8(0x0f);

forceinline __m256i popcount_avx2_64(const __m256i vec) noexcept
	{
	const __m256i lo  = _mm256_and_si256(vec, low_mask);
	const __m256i hi  = _mm256_and_si256(_mm256_srli_epi16(vec, 4), low_mask);
	const __m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
	const __m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
	const __m256i local = _mm256_add_epi8(popcnt2, popcnt1);
	const __m256i ret =_mm256_sad_epu8 (local, _mm256_setzero_si256());
	return ret;
	}

size_t at_filterAndSortByHammingDistance(uint64_t a, const std::vector<uint64_t> &b, uint64_t maxDistance, uint64_t *filteredB)
	{
	uint64_t *current_answer = filteredB;

	__m256i key = _mm256_set1_epi64x(a);
	__m256i threshold = _mm256_set1_epi64x(maxDistance);

	const uint64_t *end = &b[b.size()];
	for (uint64_t *current = (uint64_t *)&b[0]; current < end; current += 4)
		{
		__m256i data = _mm256_loadu_si256((__m256i_u *)current);
		__m256i xorResult = _mm256_xor_epi64(key, data);
		__m256i counts = popcount_avx2_64(xorResult);
//		__m256i counts = hs_popcount(xorResult);

		__mmask8 triggers = _mm256_cmple_epi64_mask(counts, threshold);
		_mm256_mask_compressstoreu_epi64(current_answer, triggers, data);
		current_answer += __popcnt16(triggers);
		}

	std::sort(filteredB, current_answer);

	return current_answer - filteredB;
	}

#endif


// Function to convert a 20mer sequence to its packed representation
uint64_t pack20mer(const char *sequence)
	{
	uint64_t packedSequence = 0;
	for (const char *end = sequence + 20; sequence < end; sequence++)
		{
		uint64_t baseEncoding;
		switch (*sequence)
			{
			case 'A':
				baseEncoding = 1;
				break;
			case 'C':
				baseEncoding = 2;
				break;
			case 'G':
				baseEncoding = 4;
				break;
			case 'T':
				baseEncoding = 7;
				break;
			default:
				baseEncoding = 0;
				break;
			}
		packedSequence = (packedSequence << 3) | baseEncoding;
		}
	return packedSequence;
	}

uint64_t *read_search_keys(const char *filename, uint64_t *key_count)
	{
	char *data = read_entire_file(filename);
	char **lines = buffer_to_list(data, key_count);

	uint64_t *keys = new uint64_t[*key_count];
	for (size_t which = 0; which < *key_count; which++)
		keys[which] = pack20mer(lines[which]);

	return keys;
	}

int main()
	{
	#ifdef __AVX512F__
		puts("AVX512");
	#else
		puts("scalar");
	#endif

	uint64_t a = 12345;
	a = 164703072086692425;

	std::vector<uint64_t> b = {54321, 67890, 98765, 24680};
	uint64_t maxDistance = 4;

	char *data = read_entire_file("OryzaPacked20mers.txt");
	uint64_t line_count;
	char **lines = buffer_to_list(data, &line_count);
	b.resize(line_count);
	for (size_t which = 0; which < line_count; which++)
		b[which] = atoll(lines[which]);

	line_count -= line_count % 4;		// so that we have a full number of SIMD registers


	uint64_t key_count;
	uint64_t *search_keys = read_search_keys("OSativaCandidates.txt", &key_count);

key_count = key_count > 1000 ? 1000 : key_count;

	auto search_time_stopwatch = JASS::timer::start();

	std::vector<uint64_t> filteredB;
	uint64_t *result_set = new uint64_t [b.size()];
	for (uint64_t key = 0; key < key_count; key++)
		{
		a = search_keys[key];
	#ifdef __AVX512F__
		auto stopwatch = JASS::timer::start();
		size_t found = at_filterAndSortByHammingDistance(a, b, maxDistance, result_set);
		auto took = JASS::timer::stop(stopwatch);
		std::cout << "Took:" << took.milliseconds() << "\n";

		// Print the filtered and sorted results
		uint64_t *end = result_set + found;
		for (uint64_t *current  = result_set; current < end; current++)
			std::cout << *current << " ";

		std::cout << std::endl;
	#else
		filteredB.clear();
		auto stopwatch = JASS::timer::start();
		filterAndSortByHammingDistance(a, b, maxDistance, filteredB);
		auto took = JASS::timer::stop(stopwatch);
		std::cout << "Took:" << took.milliseconds() << "\n";

		// Print the filtered and sorted results
		for (auto value : filteredB)
			std::cout << value << " ";

		std::cout << std::endl;
	#endif
		}

	auto took = JASS::timer::stop(search_time_stopwatch);
	std::cout << "TOTAL Took:" << took.milliseconds() << "\n";


	return 0;
	}

