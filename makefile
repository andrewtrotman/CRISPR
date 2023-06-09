all : at_hammingDist popc

at_hammingDist : at_hammingDist.cpp 
	g++ -std=c++14 -Wall -O3 -mavx2 -mavx512f -march=native at_hammingDist.cpp -o at_hammingDist

popc : popc.cpp
	g++ -std=c++14 -Wall -O3 -march=native popc.cpp -o popc


clean :
	rm popc at_hammingDist
