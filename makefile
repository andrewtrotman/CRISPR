all : at_hammingDist popc

at_hammingDist : at_hammingDist.cpp 
	g++ -O3 -mavx512f -mavx2 -march=native at_hammingDist.cpp -o at_hammingDist

popc : popc.cpp
	g++ -O3 -mavx512f -mavx2 -march=native popc.cpp -o popc


clean :
	rm popc at_hammingDist
