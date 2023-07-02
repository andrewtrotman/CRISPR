all : fast_binary_search

CPPFLAGS = -std=c++11 -pthread -O3 -Wall
CPP = g++

fast_binary_search : fast_binary_search.o score_mit_local.o
	$(CPP) $(CPPFLAGS) $^ -o $@

fast_binary_search.o : fast_binary_search.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

score_mit_local.o : score_mit_local.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@
	

clean :
	rm fast_binary_search
