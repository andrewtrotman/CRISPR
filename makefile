all : fast_binary_search extract_guides

CPPFLAGS = -std=c++17 -pthread -O3 -Wall -march=native 
CPP = g++

extract_guides : extract_guides.o file.o asserts.o
	$(CPP) $(CPPFLAGS) $^ -o $@

extract_guides.o : extract_guides.cpp file.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

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

file.o : file.cpp file.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

asserts.o : asserts.cpp asserts.h 
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean :
	rm fast_binary_search extract_guides *.o
