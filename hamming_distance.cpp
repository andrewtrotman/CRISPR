#if defined(__APPLE__) || defined(__linux__)
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>

	#define __popcnt16 __builtin_popcount
	#define __popcnt64 __builtin_popcountll
#endif

	/*
		HAMMING_DISTANCE()
		------------------
	*/
	inline uint64_t hamming_distance(uint64_t a, uint64_t b)
		{
		return __popcnt64(a ^ b);
		}

	/*
		COMPUTE_HAMMING_SET()
		---------------------
	*/
	void compute_hamming_set(uint64_t twice_the_max_distance, uint64_t key, std::vector<const uint64_t *> &positions)
		{
		const uint64_t *end = data + data_length;

		for (const uint64_t *which = data; which < end; which++)
			if (hamming_distance(key, *which) <= twice_the_max_distance)
				positions.push_back(which);
		}
