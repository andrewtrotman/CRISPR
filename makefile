all : at_hammingDist

at_hammingDist : at_hammingDist.cpp 
	g++ -std=c++14 -O3 -Wall -mavx2 -mavx512f -march=native at_hammingDist.cpp -o at_hammingDist

clean :
	rm at_hammingDist
