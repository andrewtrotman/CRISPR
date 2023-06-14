all : at_hammingDist at_getVariations getVariations

at_hammingDist : at_hammingDist.cpp 
	g++ -std=c++14 -O3 -Wall -mavx2 -mavx512f -march=native at_hammingDist.cpp -o at_hammingDist

at_getVariations : at_getVariations.cpp 
	g++ -std=c++14 -O3 -Wall -march=native at_getVariations.cpp -o at_getVariations

getVariations : getVariations.cpp 
	g++ -std=c++14 -O3 -Wall -march=native getVariations.cpp -o getVariations

clean :
	rm at_hammingDist at_getVariations getVariations
