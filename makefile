all : fast_binary_search

fast_binary_search : fast_binary_search.cpp 
	g++ -std=c++11 -pthread -O3 -Wall fast_binary_search.cpp -o fast_binary_search

clean :
	rm fast_binary_search
