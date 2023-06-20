#include <iostream>

/*
	BINARY_SEARCH()
	---------------
*/
int64_t binary_search(uint64_t *array, uint64_t elements, uint64_t key)
	{
	uint64_t lower = -1;
	uint64_t upper = elements - 1;

	while (lower + 1 != upper)
		{
		uint64_t middle = (lower + upper) / 2;
		if (array[middle] < key)
			lower = middle;
		else
			upper = middle;
		}

	/*
		The first occurance of key is at position upper if key is present in the array.
		Else, if upper > elements then we've gone off the top of the array.
	*/
	uint64_t answer = upper;
	if ((answer > elements) || (array[answer] != key))
		answer = -1;

	return answer;
	}

/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	std::cout << (1/2) << "\n";

	uint64_t array[] = {3, 4, 5, 6, 7, 8, 9};
	int n = sizeof(array) / sizeof(array[0]);

	for (int x = 0; x < 15; x++)
		std::cout << x << ":" << binary_search(array, n, x) << "\n";

	return 0;
	}

#ifdef NEVER
int main(void)
	{
	int array[] = {3, 4, 5, 6, 7, 8, 9};
	int n = sizeof(array) / sizeof(array[0]);
	int x = 4;
	int result = binarySearch(array, x, 0, n - 1);
	if (result == -1)
		printf("Not found");
	else
		printf("Element is found at index %d", result);
	return 0;
	}
#endif
