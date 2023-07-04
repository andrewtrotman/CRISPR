all : fast_binary_search

CPPFLAGS = -std=c++14 -pthread -O3 -Wall -march=native 
CPP = g++

fast_binary_search : main.o fast_binary_search.o fast_binary_search_avx512.o score_mit_local.o encode_kmer_2bit.o encode_kmer_3bit.o file.o asserts.o 
	$(CPP) $(CPPFLAGS) $^ -o $@

main.o : main.cpp finder.h hamming_distance.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

fast_binary_search.o : fast_binary_search.cpp fast_binary_search.h finder.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

fast_binary_search_avx512.o : fast_binary_search_avx512.cpp fast_binary_search_avx512.h fast_binary_search.h finder.h
	$(CPP) $(CPPFLAGS) -mavx512f -mavx2 -c $< -o $@

score_mit_local.o : score_mit_local.cpp score_mit_local.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

encode_kmer_2bit.o : encode_kmer_2bit.cpp encode_kmer_2bit.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

encode_kmer_3bit.o : encode_kmer_3bit.cpp encode_kmer_3bit.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

file.o : file.cpp file.h file_read_only.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

asserts.o : asserts.cpp asserts.h 
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean :
	rm fast_binary_search *.o
