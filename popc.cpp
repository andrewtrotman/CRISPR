#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#include <iostream>

__m256i popcount(const __m256i v)
{
  const __m256i m1 = _mm256_set1_epi8(0x55);
  const __m256i m2 = _mm256_set1_epi8(0x33);
  const __m256i m4 = _mm256_set1_epi8(0x0F);

  const __m256i t1 = _mm256_sub_epi8(v,       (_mm256_srli_epi16(v,  1) & m1));
  const __m256i t2 = _mm256_add_epi8(t1 & m2, (_mm256_srli_epi16(t1, 2) & m2));
  const __m256i t3 = _mm256_add_epi8(t2, _mm256_srli_epi16(t2, 4)) & m4;
  return _mm256_sad_epu8(t3, _mm256_setzero_si256());
}

int main(void)
{
uint64_t data[4];
uint64_t output[4];

for (int i = 0; i < sizeof(data) / sizeof(*data); i++)
	data[i] = 128 + i;

auto ans = popcount(*(__m256i *)data);
memcpy(output, &ans, sizeof(ans));

for (int i = 0; i < sizeof(output) / sizeof(*output); i++)
	printf("%02llX ", data[i]);

puts("");

for (int i = 0; i < sizeof(output) / sizeof(*output); i++)
	printf("%02llX ", output[i]);

puts("");
return 0;
}

