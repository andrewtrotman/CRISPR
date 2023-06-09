#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cstdint>

int hammingDistance(uint64_t a, uint64_t b) {
    uint64_t xorResult = a ^ b;
    return std::bitset<64>(xorResult).count();
}

void filterAndSortByHammingDistance(uint64_t a, const std::vector<uint64_t>& b, double maxDistance, std::vector<uint64_t>& filteredB) {
    // calculate the distance between a and every element in b, and return elements where dist(a,b)<maxDistance 
    filteredB.reserve(b.size());

    for (size_t i = 0; i < b.size(); i++) {
        int distance = hammingDistance(a, b[i]);
        if (distance <= maxDistance) {
            filteredB.push_back(b[i]);
        }
    }

    std::sort(filteredB.begin(), filteredB.end());
}

int main() {
    uint64_t a = 12345;
    std::vector<uint64_t> b = {54321, 67890, 98765, 24680};
    double maxDistance = 2.0;

    std::vector<uint64_t> filteredB;
    filterAndSortByHammingDistance(a, b, maxDistance, filteredB);

    // Print the filtered and sorted results
    for (const auto& value : filteredB) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}

