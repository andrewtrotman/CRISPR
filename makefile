all : at_hammingDist at_getVariations getVariations at_evalGuides fast_binary_search avx512_binary_search

at_hammingDist : at_hammingDist.cpp 
	g++ -std=c++14 -pthread -O3 -Wall -mavx2 -mavx512f -march=native at_hammingDist.cpp -o at_hammingDist

fast_binary_search : fast_binary_search.cpp 
	g++ -std=c++11 -pthread -O3 -Wall fast_binary_search.cpp -o fast_binary_search

avx512_binary_search: avx512_binary_search.cpp 
	g++ -std=c++14 -pthread -O3 -Wall -mavx2 -mavx512f -march=native avx512_binary_search.cpp -o avx512_binary_search

at_getVariations : at_getVariations.cpp 
	g++ -std=c++17 -O3 -Wall -mavx2 -mavx512f -march=native at_getVariations.cpp -o at_getVariations

getVariations : getVariations.cpp 
	g++ -std=c++17 -O3 -Wall -mavx2 -mavx512f -march=native getVariations.cpp -o getVariations

at_evalGuides : at_evalGuides.cpp
	g++ -std=c++17 -O3 -Wall -march=native at_evalGuides.cpp -o at_evalGuides


clean :
	rm at_hammingDist at_getVariations getVariations at_evalGuides fast_binary_search avx512_binary_search
