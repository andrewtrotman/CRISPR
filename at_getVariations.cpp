//
//  main.cpp
//  twentyMerVariations
//
//  Created by Shlomo Geva on 9/6/2023.
//

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <new>
#include <immintrin.h>


// Function to convert a 20mer sequence to its packed representation
uint64_t pack20mer(const std::string& sequence) {
    uint64_t packedSequence = 0;
    for (int i = 0; i < 20; i++) {
        uint64_t baseEncoding;
        switch (sequence[i]) {
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

// Function to convert a packed uint64_t sequence to its 20mer representation
std::string unpack20mer(uint64_t packedSequence) {
    std::string sequence;
    for (int i = 19; i >= 0; i--) {
        uint64_t baseEncoding = (packedSequence >> (i * 3)) & 7;
        char base;
        switch (baseEncoding) {
            case 1:
                base = 'A';
                break;
            case 2:
                base = 'C';
                break;
            case 4:
                base = 'G';
                break;
            case 7:
                base = 'T';
                break;
            default:
                base = 'N';
                break;
        }
        sequence += base;
    }
    return sequence;
}

std::vector<std::vector<int>> generatePositionSets() {
    std::vector<std::vector<int>> positionSets;
    for (int i = 1; i <= 17; i++) {
        for (int j = i + 1; j <= 18; j++) {
            for (int k = j + 1; k <= 19; k++) {
                for (int l = k + 1; l <= 20; l++) {
                    positionSets.push_back({i, j, k, l});
                }
            }
        }
    }
    return positionSets;
}

std::vector<std::string> generateChoicesWithReplacement() {
    std::vector<std::string> choices;
    std::string bases = "ACGT";

    for (char c1 : bases) {
        for (char c2 : bases) {
            for (char c3 : bases) {
                for (char c4 : bases) {
                    std::string choice(4, ' ');
                    choice[0] = c1;
                    choice[1] = c2;
                    choice[2] = c3;
                    choice[3] = c4;
                    choices.push_back(choice);
                }
            }
        }
    }

    return choices;
}

std::vector<std::string> generateMasks(const std::vector<std::vector<int>>& positionSets, const std::vector<std::string>& choices) {
    std::vector<std::string> masks;

    for (const auto& positionSet : positionSets) {
        for (const auto& choice : choices) {
            std::string mask(20, '.');

            for (unsigned int j = 0; j < positionSet.size(); ++j) {
                int position = positionSet[j];
                char base = choice[j];
                mask[position - 1] = base;
            }

            masks.push_back(mask);
        }
    }

    return masks;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> packMasks(const std::vector<std::string>& masks) {
    std::vector<uint64_t> changeMasks;
    std::vector<uint64_t> preserveMasks;

    for (const auto& mask : masks) {
        uint64_t changeMask = 0;
        uint64_t preserveMask = 0x0FFFFFFFFFFFFFFFULL;

        for (unsigned int i = 0; i < mask.size(); ++i) {
            char base = mask[i];

            if (base != '.') {
                int position = (19 - i) * 3;
                switch (base) {
                    case 'A':
                        changeMask |= (1ULL << position);
                        preserveMask &= ~(7ULL << position);
                        break;
                    case 'C':
                        changeMask |= (2ULL << position);
                        preserveMask &= ~(7ULL << position);
                        break;
                    case 'G':
                        changeMask |= (4ULL << position);
                        preserveMask &= ~(7ULL << position);
                        break;
                    case 'T':
                        changeMask |= (7ULL << position);
                        preserveMask &= ~(7ULL << position);
                        break;
                }
            }
        }

        changeMasks.push_back(changeMask);
        preserveMasks.push_back(preserveMask);
    }

    return {changeMasks, preserveMasks};
}

size_t getVariations64(uint64_t *variations, const std::vector<uint64_t>& changeMasks, const std::vector<uint64_t>& preserveMasks, const uint64_t kmer)
	{
	uint64_t *into = variations;
	for (std::size_t i = 0; i < changeMasks.size(); ++i)
		{
		const uint64_t variation = (kmer & preserveMasks[i]) | changeMasks[i];
		*into++ = variation;
		}

	// Sort and remove duplicates
	std::sort(variations, into);
	into = std::unique(variations, into);

	return into - variations;
	}

size_t getVariations256(uint64_t *variations, const std::vector<uint64_t>& changeMasks, const std::vector<uint64_t>& preserveMasks, const uint64_t kmer)
	{
	__m256i *into = (__m256i *)variations;
	__m256i kmer_256 = _mm256_set1_epi64x(kmer);
	__m256i *cm, *pm;

	__m256i *end = (__m256i *)&preserveMasks[preserveMasks.size()];
	for (pm = (__m256i *)&preserveMasks[0], cm = (__m256i *)&changeMasks[0]; pm < end; pm++, cm++)
		{
		__m256i preserve_mask = _mm256_load_si256(pm);
		__m256i change_mask = _mm256_load_si256(cm);

		__m256i variation = _mm256_or_si256(_mm256_and_si256(kmer_256, preserve_mask), change_mask);
		_mm256_store_si256(into, variation);

		into++;
		}

	// Sort and remove duplicates
	std::sort(variations, (uint64_t *)into);
	return std::unique(variations, (uint64_t *)into) - variations;
	}


size_t getVariations512(uint64_t *variations, const std::vector<uint64_t>& changeMasks, const std::vector<uint64_t>& preserveMasks, const uint64_t kmer)
	{
	__m512i *into = (__m512i *)variations;
	__m512i kmer_512 = _mm512_set1_epi64(kmer);
	__m512i *cm, *pm;

	__m512i *end = (__m512i *)&preserveMasks[preserveMasks.size()];
	for (pm = (__m512i *)&preserveMasks[0], cm = (__m512i *)&changeMasks[0]; pm < end; pm++, cm++)
		{
		__m512i preserve_mask = _mm512_load_si512(pm);
		__m512i change_mask = _mm512_load_si512(cm);

		__m512i variation = _mm512_or_si512(_mm512_and_si512(kmer_512, preserve_mask), change_mask);
		_mm512_store_si512(into, variation);

		into++;
		}

	// Sort and remove duplicates
	std::sort(variations, (uint64_t *)into);
	return std::unique(variations, (uint64_t *)into) - variations;
	}

int usage(const char *exename)
	{
	printf("Usage:%s <64 | 256 | 512>", exename);
	return 1;
	}

int main(int argc, const char *argv[])
	{
	if (argc != 2)
		exit(usage(argv[0]));

	long bits_wide;
	bits_wide = atol(argv[1]);
	if (bits_wide != 64 && bits_wide != 256 && bits_wide != 512)
		exit(usage(argv[0]));

	std::cout << "Generating masks" << std::endl;
	auto start = std::chrono::high_resolution_clock::now(); // Start timing

	std::vector<std::vector<int>> positionSets = generatePositionSets();
	std::vector<std::string> choices = generateChoicesWithReplacement();

	std::vector<std::string> masks = generateMasks(positionSets, choices);
	std::pair<std::vector<uint64_t>, std::vector<uint64_t>> packedMasks = packMasks(masks);
	std::vector<uint64_t> changeMasks = packedMasks.first;
	std::vector<uint64_t> preserveMasks = packedMasks.second;

	auto end = std::chrono::high_resolution_clock::now(); // End timing
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken: " << duration << " microseconds" << std::endl;

	std::cout << "Number of masks: " << masks.size() << std::endl;
	std::cout << "Number of packed masks: " << changeMasks.size() << std::endl;

	std::string kmer = "ACGTTGCATTAAGGCCGGAC";
	uint64_t packedKmer = pack20mer(kmer);

	start = std::chrono::high_resolution_clock::now(); // Start timing

	auto method = getVariations64;
	if (bits_wide == 64)
		method = getVariations64;
	else if (bits_wide == 256)
		method = getVariations256;
	else if (bits_wide == 512)
		method = getVariations512;

	__m512i *variations_aligned = new __m512i[masks.size() / 8 + 1];
	uint64_t *variations = (uint64_t *)variations_aligned;


	size_t variations_size = method(variations, changeMasks, preserveMasks, packedKmer);
	
	end = std::chrono::high_resolution_clock::now(); // End timing
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Number of variations: " << variations_size << std::endl;
	for (int i=0; i<10; i++)
		std::cout << unpack20mer(variations[i]) << std::endl;

	std::cout << "Time taken: " << duration << " microseconds" << std::endl;

	delete [] variations;
	return 0;
	}
