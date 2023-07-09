all : fast_binary_search extract_guides

CPPFLAGS = -std=c++17 -pthread -O3 -Wall -march=native 
CXX = g++

extract_guides : extract_guides.o file.o asserts.o
	$(CXX) $(CPPFLAGS) $^ -o $@

fast_binary_search : main.o fast_binary_search.o score_mit_local.o encode_kmer_2bit.o encode_kmer_3bit.o file.o asserts.o
	$(CXX) $(CPPFLAGS) $^ -o $@

asserts.o: asserts.cpp asserts.h

encode_kmer_2bit.o: encode_kmer_2bit.cpp encode_kmer_2bit.h

encode_kmer_3bit.o: encode_kmer_3bit.cpp encode_kmer_3bit.h

extract_guides.o: extract_guides.cpp file.h

fast_binary_search.o: fast_binary_search.cpp forceinline.h score_mit_local.h fast_binary_search.h job.h finder.h sequence_score_pair.h encode_kmer_2bit.h

fast_binary_search_avx512.o: fast_binary_search_avx512.cpp fast_binary_search_avx512.h fast_binary_search.h job.h finder.h sequence_score_pair.h score_mit_local.h encode_kmer_2bit.h

file.o: file.cpp file.h asserts.h

main.o: main.cpp file.h finder.h job.h sequence_score_pair.h encode_kmer_2bit.h encode_kmer_3bit.h fast_binary_search.h score_mit_local.h

score_mit_local.o: score_mit_local.cpp score_mit_local.h



clean :
	rm fast_binary_search extract_guides *.o
