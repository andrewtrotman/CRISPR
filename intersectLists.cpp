#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>

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

void compute_intersection_list(const uint64_t* A, size_t sizeA, const uint64_t* B, size_t sizeB, std::vector<uint64_t>& matches, std::vector<size_t>& positions) {
    // intersect arrays A and B, returning matching elements in 'matches'
    // and the respective match  positions in 'position'
    //
    size_t i = 0;
    size_t j = 0;

    while (i < sizeA && j < sizeB) {
        if (A[i] == B[j]) {
            matches.push_back(A[i]);
            positions.push_back(j + 1);  // Add 1 to match MATLAB's 1-based indexing
            i++;
            j++;
        } else if (A[i] < B[j]) {
            i = bisect_left(A, sizeA, B[j], i + 1);
        } else {
            j = bisect_left(B, sizeB, A[i], j + 1);
        }
    }
}

int main() {
    // Sample data
    const uint64_t A[] = {1, 2, 3, 4, 5};
    size_t sizeA = sizeof(A) / sizeof(A[0]);
    const uint64_t B[] = {2, 4, 6, 8};
    size_t sizeB = sizeof(B) / sizeof(B[0]);

    std::vector<uint64_t> matches;
    std::vector<size_t> positions;

    compute_intersection_list(A, sizeA, B, sizeB, matches, positions);

    // Print the results
    std::cout << "Matches: ";
    for (const auto& match : matches) {
        std::cout << match << " ";
    }
    std::cout << std::endl;

    std::cout << "Positions: ";
    for (const auto& position : positions) {
        std::cout << position << " ";
    }
    std::cout << std::endl;

    return 0;
}

