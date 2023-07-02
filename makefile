all : fast_binary_search

CPPFLAGS = -std=c++11 -pthread -O3 -Wall
CPP = g++

fast_binary_search : main.o fast_binary_search.o score_mit_local.o encode_kmer_2bit.o encode_kmer_3bit.o file.o file_read_only.o
	$(CPP) $(CPPFLAGS) $^ -o $@

main.o : main.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

fast_binary_search.o : fast_binary_search.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

score_mit_local.o : score_mit_local.cpp score_mit_local.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

encode_kmer_2bit.o : encode_kmer_2bit.cpp encode_kmer_2bit.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

encode_kmer_3bit.o : encode_kmer_3bit.cpp encode_kmer_3bit.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

file.o : file.cpp file.h file_read_only.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

file_read_only.o : file_read_only.cpp file_read_only.h
	$(CPP) $(CPPFLAGS) -c $< -o $@

clean :
	rm fast_binary_search *.o
