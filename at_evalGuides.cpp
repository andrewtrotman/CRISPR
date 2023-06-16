#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <string>
#include <sstream>
#include <random>

// Declare mutexes for synchronization
std::mutex matchesMutex;
std::mutex positionsMutex;

uint64_t encodingTable[256];

// Function to perform intersection between two uint64_t arrays
template <typename T>
size_t bisect_left(const T* arr, size_t size, const T& value, size_t start = 0) {
    size_t lo = start;
    size_t hi = size;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (arr[mid] < value) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// The bisect_right C++ function is used to perform binary search on a sorted sequence
//   and find the index at which a given value should be inserted while maintaining
//   the sorted order.  Working in the opposite direction to the library function bisect_left

template <typename T>
size_t bisect_right(const T* arr, size_t size, const T& value, size_t start) {
    size_t left = start;
    size_t right = size;

    while (left < right) {
        size_t mid = left + (right - left) / 2;

        if (value < arr[mid]) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

/* *** KEEP ***
// Function to perform intersection between two uint64_t arrays
void compute_intersection_list(const uint64_t* A, size_t sizeA, const uint64_t* B, size_t sizeB, std::vector<uint64_t>& matches, std::vector<size_t>& positions) {
    size_t i = 0;
    size_t j = 0;

    while (i < sizeA && j < sizeB) {
        if (A[i] == B[j]) {
            matches.push_back(A[i]);
            positions.push_back(j + 1); // Add 1 to match vi editor 1-based line indexing
            i++;
            j++;
        } else if (A[i] < B[j]) {
            i = bisect_left(A, sizeA, B[j], i + 1);
        } else {
            j = bisect_left(B, sizeB, A[i], j + 1);
        }
    }
}
*/

///*
// Function to perform intersection between two uint64_t arrays
void compute_intersection_list(const uint64_t* A, size_t sizeA, const uint64_t* B, size_t sizeB, std::vector<uint64_t>& matches, std::vector<size_t>& positions) {
    size_t startA = 0;
    size_t startB = 0;
    size_t endA = sizeA;
    size_t endB = sizeB;

   while (true) {
       // lookup from the left
       if (A[startA] == B[startB]) {
           matches.push_back(A[startA]);
           positions.push_back(startB + 1); // Add 1 to match vi editor 1-based line indexing (just convenient)
           startA++;
           startB++;
       } else if (A[startA] < B[startB]) {
           startA = bisect_left(A, sizeA, B[startB], startA + 1);
       } else {
           startB = bisect_left(B, sizeB, A[startA], startB + 1);
       }
       if (!(startA < sizeA && startB < sizeB))
           break;

       // lookup from the right
       if (A[endA] == B[endB]) {
           matches.push_back(A[endA]);
           positions.push_back(endB + 1); // Add 1 to match vi editor 1-based line indexing (just convenient)
           endA--;
           endB--;
       } else if (A[endA] > B[endB]) {
           endA = bisect_right(A, sizeA, B[endB], endA - 1);
       } else {
           endB = bisect_left(B, sizeB, A[endA], endB + 1);
       }
       if (!(endA >= 0 && endB >= 0))
           break;
       

   }
}
/*
void compute_intersection_list(const uint64_t* A, size_t sizeA, const uint64_t* B, size_t sizeB, std::vector<uint64_t>& matches, std::vector<size_t>& positions) {
    size_t i = 0;
    size_t j = 0;

    while (i < sizeA && j < sizeB) {
        if (A[i] == B[j]) {
            matches.push_back(A[i]);
            positions.push_back(j + 1);  // Add 1 to match vi's 1-based line indexing
            i++;
            j++;
        } else if (A[i] < B[j]) {
            i = bisect_left(A, sizeA, B[j], i + 1);
        } else {
            j = bisect_left(B, sizeB, A[i], j + 1);
        }
    }
}
*/

//// Function to convert a 20mer sequence to its packed representation
//uint64_t pack20mer(const std::string& sequence) {
//    uint64_t packedSequence = 0;
//    for (int i = 0; i < 20; i++) {
//        uint64_t baseEncoding;
//        switch (sequence[i]) {
//            case 'A':
//                baseEncoding = 1;
//                break;
//            case 'C':
//                baseEncoding = 2;
//                break;
//            case 'G':
//                baseEncoding = 4;
//                break;
//            case 'T':
//                baseEncoding = 7;
//                break;
//            default:
//                baseEncoding = 0;
//                break;
//        }
//        packedSequence = (packedSequence << 3) | baseEncoding;
//    }
//    return packedSequence;
//}

//uint64_t pack20mer(const std::string& sequence) {
//    static const uint64_t encodingTable[256] = {
//        // Initialize the lookup table with base encodings
//        // 0 for invalid bases, 1 for 'A', 2 for 'C', 4 for 'G', 7 for 'T'
//        // Note: Make sure to initialize the unused entries to 0
//        // for (int i = 0; i < 256; i++) encodingTable[i] = 0; (if necessary)
//        ['A'] = 1,
//        ['C'] = 2,
//        ['G'] = 4,
//        ['T'] = 7,
//    };
//
//    uint64_t packedSequence = 0;
//    for (int i = 0; i < 20; i++) {
//        uint64_t baseEncoding = encodingTable[sequence[i]];
//        packedSequence = (packedSequence << 3) | baseEncoding;
//    }
//    return packedSequence;
//}

uint64_t pack20mer(const std::string& sequence) {
    uint64_t packedSequence = 0;
    for (int i = 0; i < 20; i++) {
        uint64_t baseEncoding = encodingTable[sequence[i]];
        packedSequence = (packedSequence << 3) | baseEncoding;
    }
    return packedSequence;
}

uint64_t pack20mer(const char *sequence)
	{
	uint64_t packedSequence = 0;
	for (int i = 0; i < 20; i++)
		{
		uint64_t baseEncoding = encodingTable[sequence[i]];
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
// Function to load sequences from a file into a vector
#ifdef NEVER
	std::vector<std::string> loadSequencesFromFile(const std::string& filename) {
		 std::vector<std::string> sequences;
		 std::ifstream file(filename);
		 if (file.is_open()) {
			  std::string line;
			  while (std::getline(file, line)) {
					if (line.size() == 20) {
						 sequences.push_back(line);
					}
			  }
			  file.close();
		 }
		 return sequences;
	}
#else
	/*
		LOADSEQUENCESFROMFILE()
		-----------------------
	*/
	std::vector<std::string> loadSequencesFromFile(const std::string& filename)
		{
		std::vector<std::string> sequences;
		char *data = read_entire_file(filename.c_str());
		char *start = data;
		char *eoln;
		if (data != NULL)
			{
			eoln = strchr(start, '\n');
			if (eoln != NULL)
				if (eoln - start == 20)
					sequences.push_back(std::string(start, 20));
			}
		free(data);
		return sequences;
		}
#endif


std::vector<uint64_t> loadPackedGenomeGuidesFromFile(const std::string& filename) {
    std::vector<uint64_t> packedGenomeGuides;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        uint64_t packedGuide = std::stoull(line);
        packedGenomeGuides.push_back(packedGuide);
    }
    return packedGenomeGuides;
}

#ifndef NEVER
	std::vector<uint64_t> loadGuidesFromFile(const std::string& filename) {
		 std::vector<uint64_t> packedGenomeGuides;
		 std::ifstream guideFile(filename);
		 if (!guideFile.is_open()) {
			  std::cerr << "Error opening guide file." << std::endl;
			  exit(1);
		 }
    
		 std::string line;
		 while (std::getline(guideFile, line)) {
			  std::istringstream iss(line);
			  std::string guide;
			  uint64_t packedGuide;
			  int sequenceNumber, startPosition;
	//        std::cout << line << std::endl;
        
			  if (iss >> guide >> sequenceNumber >> startPosition) {
					// Process the extracted guide, sequenceNumber, and startPosition
					packedGuide = pack20mer(guide);
					packedGenomeGuides.push_back(packedGuide);
			  } else {
					std::cerr << "Error parsing line: " << line << std::endl;
					exit(2);
			  }
		  }
		 std::sort(packedGenomeGuides.begin(),packedGenomeGuides.end());
		 return packedGenomeGuides;
	}
#else
	std::vector<uint64_t> loadGuidesFromFile(const std::string& filename)
		{
		std::vector<uint64_t> packedGenomeGuides;
		char *guide;
		char *space;
		char *data = read_entire_file(filename.c_str());
		if (data == NULL)
			{
			std::cerr << "Error opening guide file." << std::endl;
			exit(1);
			}
		
		guide = data - 1;
		do
			{
			guide++;
			packedGenomeGuides.push_back(pack20mer(guide));
			guide = strchr(guide, '\n');
			}
		while (guide != NULL && *(guide + 1) != '\0');

		std::sort(packedGenomeGuides.begin(),packedGenomeGuides.end());
		free(data);
		return packedGenomeGuides;
		}

#endif

#ifdef NEVER
// Function to generate variations with 0 to 4 replacements away from a given 20mer
void generateVariations(std::string& sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0) {
    if (replacements > 4) {
        return;
    }
    char originalBase;
    variations.push_back(pack20mer(sequence));
    for (int i = position; i < 20; i++) {
        for (char base : {'A', 'C', 'G', 'T'}) {
            if (sequence[i] != base) {
                originalBase = sequence[i];  // Store the original base value
                sequence[i] = base;  // Modify the base in 'modifiedSequence'
                generateVariations(sequence, variations, replacements + 1, i + 1);
                sequence[i] = originalBase;  // Restore the original base value
            }
        }
    }
}
#else
// Function to generate variations with 0 to 4 replacements away from a given 20mer
void generateVariations_binary(uint64_t sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0)
	{
	if (replacements > 4)
		return;

	variations.push_back(sequence);
	for (uint64_t i = position; i < 20; i++)
		{
		uint64_t was = (sequence >> (i * 3)) & 7;
		uint64_t knock_out = (sequence & ~(7ULL << (i * 3)));
		for (uint64_t base : {1ULL, 2ULL, 4ULL, 7ULL})
			if (was != base)
				{
				uint64_t new_sequence = knock_out | (base << (i * 3));
				generateVariations_binary(new_sequence, variations, replacements + 1, i + 1);
				}
		}
	}
void generateVariations(std::string& sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0)
	{
	generateVariations_binary(pack20mer(sequence), variations, replacements, position);
	}



	#ifdef NEVER
	// Function to generate variations with 0 to 4 replacements away from a given 20mer
	void generateVariations(std::string& sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0)
	{
	puts("ASPT generateVariations");
	std::vector<std::vector<int>> positionSets;

	uint64_t encoded_sequence = pack20mer(sequence);
	uint64_t new_sequence;
	const uint64_t base_encodings[] = {1, 2, 4, 7};//        // 0 for invalid bases, 1 for 'A', 2 for 'C', 4 for 'G', 7 for 'T'

	for (int i = 0; i < 17; i++)
		{
		uint64_t new_mask = encoded_sequence & ~(7ULL << i * 3);
		for (int base = 0; base < 4; base++)
			{
			uint64_t encoded_sequence = new_mask | (base_encodings[base] << i * 3);
			for (int j = i + 1; j < 18; j++)
				{
				uint64_t new_mask = encoded_sequence & ~(7ULL << j * 3);
				for (int base = 0; base < 4; base++)
					{
					uint64_t encoded_sequence = new_mask | (base_encodings[base] << j * 3);
					for (int k = j + 1; k < 19; k++)
						{
						uint64_t new_mask = encoded_sequence & ~(7ULL << k * 3);
						for (int base = 0; base < 4; base++)
							{
							uint64_t encoded_sequence = new_mask | (base_encodings[base] << k * 3);
							for (int l = k + 1; l < 20; l++)
								{
								uint64_t new_mask = encoded_sequence & ~(7ULL << l * 3);
								for (int base = 0; base < 4; base++)
									{
									uint64_t encoded_sequence = new_mask | (base_encodings[base] << l * 3);
									variations.push_back(encoded_sequence);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	#endif

#endif

// Function to generate variations with 0 to 4 replacements away from a given 20mer
int generateVariations2(std::string& sequence, std::vector<uint64_t>& variations, int index, int replacements = 0, int position = 0) {
    if (replacements > 4) {
        return index;
    }
    char originalBase;
    variations[index++]=(pack20mer(sequence));
    for (int i = position; i < 20; i++) {
        for (char base : {'A', 'C', 'G', 'T'}) {
            if (sequence[i] != base) {
                originalBase = sequence[i];  // Store the original base value
                sequence[i] = base;  // Modify the base in 'modifiedSequence'
                index = generateVariations2(sequence, variations, index,replacements + 1, i + 1);
                sequence[i] = originalBase;  // Restore the original base value
            }
        }
    }
    return index;
}

// function to select, at random, unpacked guides from the genome set of all packed guides
std::vector<std::string> selectRandomVectors(const std::vector<uint64_t>& guides, int seed, int count) {
    std::mt19937 rng(seed);  // Initialize the random number generator with the given seed
    std::vector<std::string> selected;

    if (count > guides.size()) {
        std::cout << "Warning: count is larger than vector size, returning all" << std::endl;
        count = (int) guides.size();  // Return all vectors if count exceeds or equals the size of the input vector
    }

    std::uniform_int_distribution<int> distribution(0, (int) guides.size() - 1);

    for (int i = 0; i < count; ++i) {
        int index = distribution(rng);  // Generate a random index within the range of the input vector
        selected.push_back(unpack20mer(guides[index]));
    }

    return selected;
}


int main() {
    auto startMain = std::chrono::steady_clock::now(); // Start timing    
    // Initialize global lookup table with base encodings
    // 0 for invalid bases, 1 for 'A', 2 for 'C', 4 for 'G', 7 for 'T'
    // Note: Make sure to initialize the unused entries to 0
    // for (int i = 0; i < 256; i++) encodingTable[i] = 0; (if necessary)
    encodingTable['A'] = 1;
    encodingTable['C'] = 2;
    encodingTable['G'] = 4;
    encodingTable['T'] = 7;

    std::cout << "Loading test guides and genome guides" <<  std::endl;
    /* load genome guide */
    
    //    std::string genomeGuidesFilename = "/Users/geva/Crispr/OryzaPacked20mers.txt";
    //    std::vector<uint64_t> packedGenomeGuides = loadPackedGenomeGuidesFromFile(genomeGuidesFilename);
    
    std::string GuidesFilename = "OryzaSativaGuides.txt";

    auto load_guided_from_file_start = std::chrono::steady_clock::now();
    std::vector<uint64_t> packedGenomeGuides = loadGuidesFromFile(GuidesFilename);
    auto load_guided_from_file_end = std::chrono::steady_clock::now();
    auto load_guided_from_file_duration = std::chrono::duration_cast<std::chrono::microseconds>(load_guided_from_file_end - load_guided_from_file_start).count();
    std::cout << "loadGuidesFromFile: " << load_guided_from_file_duration/1000000.001 << " seconds" << std::endl;

    /* load test guides */
    
    //    std::string inputSequencesFilename = "/Users/geva/Crispr/testGuides.txt";
//    std::string inputSequencesFilename = "/Users/geva/Crispr/OSativaCandidates.txt";
    
//    std::vector<std::string> testGuides = loadSequencesFromFile(inputSequencesFilename);
//    count = (int) testGuides.size();

    // select test guides at random, from the genome guides
    int seed = 13; // for repeatability
    int TESTSIZE = 10000; //packedGenomeGuides.size(); //10000; // number of random guides to test
    std::vector<std::string> testGuides = selectRandomVectors(packedGenomeGuides, seed, TESTSIZE);
    
    std::cout << "Loaded " << testGuides.size() << " test guides, " <<  packedGenomeGuides.size() << " genome guides" << std::endl;
    
    ////////////////    /Users/geva/Crispr/humanGenome/GRCh38_latest_genomic.fna
    ///

    std::vector<std::vector<uint64_t>> allMatches;
    std::vector<std::vector<size_t>> allPositions;

    auto start2 = std::chrono::steady_clock::now(); // Start timing
    
    bool Parallelise = true;
    
    if (Parallelise)
    {
        // Define the number of threads to use for parallel execution
        int numThreads = std::thread::hardware_concurrency();
//        numThreads=128;
        numThreads=1;
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        
        // Split the inputSequences into equal-sized chunks for parallel processing
        const size_t chunkSize = testGuides.size() / numThreads;
        size_t startIndex = 0;
        
        // Function to be executed in parallel by each thread
        auto processChunk = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {

//                std::vector<uint64_t> variations(424996);
//                generateVariations2(testGuides[i], variations, 0, 0, 0);
                
                std::vector<uint64_t> variations;
				    auto generateVariations_start = std::chrono::steady_clock::now();
                generateVariations(testGuides[i], variations, 0, 0);
					 auto generateVariations_end = std::chrono::steady_clock::now();
					 auto generateVariations_duration = std::chrono::duration_cast<std::chrono::microseconds>(generateVariations_end - generateVariations_start).count();
					 std::cout << "generateVariations: " << generateVariations_duration/1000000.001 << " seconds" << std::endl;



std::cout << "BASED ON:" << testGuides[i] << "\n";
for (const auto x : variations)
	{
	std::cout << x << '\n';
	}
exit(0);


                variations.erase(variations.begin());
                std::sort(variations.begin(), variations.end());
                
                std::vector<uint64_t> matches;
                std::vector<size_t> positions;
                //            auto start4 = std::chrono::steady_clock::now(); // Start timing
                compute_intersection_list(variations.data(), variations.size(), packedGenomeGuides.data(), packedGenomeGuides.size(), matches, positions);
                //            auto end4 = std::chrono::steady_clock::now(); // End timing
                //            auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count();
                //            std::cout << "Intersection Execution time: " << duration4/1000000.001 << " seconds" << std::endl;
                
                // Lock the mutexes before modifying the shared vectors
                std::lock_guard<std::mutex> matchesLock(matchesMutex);
                std::lock_guard<std::mutex> positionsLock(positionsMutex);
                
                allMatches.push_back(matches);
                allPositions.push_back(positions);
                
                // Mutexes are automatically released when the lock_guard goes out of scope
                
            }
        };
        std::cout << "Launching " << numThreads << " threads" << std::endl;
        // Create and launch threads to process each chunk
        for (int i = 0; i < numThreads - 1; ++i) {
            threads.emplace_back(processChunk, startIndex, startIndex + chunkSize);
            startIndex += chunkSize;
        }
        // The last thread may need to process the remaining items if the chunk size is not perfectly divisible
        threads.emplace_back(processChunk, startIndex, testGuides.size());
        
        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }
    }
    else {
        // Single theared
        int variationTime=0, IntersectTime=0;
        std::cout << "Single threaded execution\n";
        for (int i=0; i<testGuides.size(); i++) {
            std::vector<uint64_t> variations;
            std::vector<uint64_t> matches;
            std::vector<size_t> positions;
            auto start = std::chrono::steady_clock::now(); // Start timing
            generateVariations(testGuides[i], variations, 0, 0);
            variations.erase(variations.begin()); // self = zero variatios
            auto end = std::chrono::steady_clock::now(); // End timing
            variationTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::sort(variations.begin(), variations.end());
            start = std::chrono::steady_clock::now(); // Start timing
            compute_intersection_list(variations.data(), variations.size(), packedGenomeGuides.data(), packedGenomeGuides.size(), matches, positions);
            end = std::chrono::steady_clock::now(); // End timing
            IntersectTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            allMatches.push_back(matches);
            allPositions.push_back(positions);
        }
        std::cout << "Variation execution time: " << variationTime/1000000/(TESTSIZE+0.0001) << " seconds" << std::endl;
        std::cout << "Intersection execution time: " << IntersectTime/1000000/(TESTSIZE+0.0001) << " seconds" << std::endl;
    }
    auto end2 = std::chrono::steady_clock::now(); // End timing
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count();

        // Calculate statistics
        size_t totalCount = allMatches.size();
        size_t emptyCount = std::count_if(allPositions.begin(), allPositions.end(), [](const std::vector<size_t>& positions) {
            return positions.empty();
        });
        double meanMatches = 0.0;
        double sumMatches = 0.0;
        size_t minMatches = std::numeric_limits<size_t>::max();
        size_t maxMatches = 0;
    
        for (const auto& matches : allMatches) {
            size_t count = matches.size();
            sumMatches += count;
            minMatches = std::min(minMatches, count);
            maxMatches = std::max(maxMatches, count);
        }
    
        meanMatches = sumMatches / totalCount;
    
        // Calculate standard deviation
        double variance = 0.0;
        for (const auto& matches : allMatches) {
            double diff = matches.size() - meanMatches;
            variance += (diff * diff);
        }
        variance /= totalCount;
        double stdDevMatches = std::sqrt(variance);
    
        // Print the results
        std::cout << "Mean number of matches: " << meanMatches << std::endl;
        std::cout << "Standard deviation of matches: " << stdDevMatches << std::endl;
        std::cout << "Minimum matches: " << minMatches << std::endl;
        std::cout << "Maximum matches: " << maxMatches << std::endl;
        std::cout << "Number of unique test guides: " << emptyCount << std::endl;
        std::cout << "Number of test guides with off-targets: " << TESTSIZE-emptyCount << std::endl;
        std::cout << "Mean Execution time (getvar+intersect): " << duration2/1000000.001 << " seconds" << std::endl;
        auto endMain = std::chrono::steady_clock::now(); // End timing
        auto durationMain = std::chrono::duration_cast<std::chrono::microseconds>(endMain - startMain).count();
       std::cout << "Total Execution time: " << durationMain/1000000.001 << " seconds" << std::endl;
    return 0;
}


