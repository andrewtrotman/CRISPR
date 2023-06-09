#include <algorithm>
#include "mex.h"
#include <iostream>
#include <vector>

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

// Function to generate variations with 0 to 4 replacements away from a given 20mer
void generateVariations(std::string& sequence, std::vector<uint64_t>& variations, int replacements = 0, int position = 0) {
    if (replacements > 4) {
        return;
    }
    variations.push_back(pack20mer(sequence));
    for (int i = position; i < 20; i++) {
        for (char base : {'A', 'C', 'G', 'T'}) {
            if (sequence[i] != base) {
                char originalBase = sequence[i];  // Store the original base value
                sequence[i] = base;  // Modify the base in 'sequence'
                generateVariations(sequence, variations, replacements + 1, i + 1);
                sequence[i] = originalBase;  // Restore the original base value
            }
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("packCRISPR20mers:invalidNumInputs", "One input argument is required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("packCRISPR20mers:invalidNumOutputs", "Too many output arguments.");
    }
    if (!mxIsChar(prhs[0]) || mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("packCRISPR20mers:invalidInput", "Input must be a single row character array.");
    }

    char* sequence = mxArrayToString(prhs[0]);
    std::string inputSequence(sequence);
    mxFree(sequence);

    std::vector<uint64_t> variations;
    generateVariations(inputSequence, variations);
    variations.erase(variations.begin());
    std::sort(variations.begin(), variations.end());

    plhs[0] = mxCreateNumericMatrix(variations.size(), 1, mxUINT64_CLASS, mxREAL);
    uint64_t* packedKmers = reinterpret_cast<uint64_t*>(mxGetData(plhs[0]));
    std::copy(variations.begin(), variations.end(), packedKmers);
}

