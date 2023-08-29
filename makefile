all : fast_binary_search.exe extract_guides.exe

CPP_LINK_FLAGS = /std:c++11 /EHsc /Zi
CPPFLAGS = $(CPP_LINK_FLAGS) /c
CPP = cl

fast_binary_search.exe : fast_binary_search.obj main.obj score_mit_local.obj encode_kmer_2bit.obj encode_kmer_3bit.obj file.obj asserts.obj
	$(CPP) $(CPP_LINK_FLAGS) $**

extract_guides.exe : extract_guides.obj file.obj asserts.obj
	$(CPP) $(CPP_LINK_FLAGS) $**

main.obj : main.cpp finder.h hamming_distance.h

fast_binary_search.obj : fast_binary_search.cpp fast_binary_search.h finder.h

score_mit_local.obj : score_mit_local.cpp score_mit_local.h

encode_kmer_2bit.obj : encode_kmer_2bit.cpp encode_kmer_2bit.h

encode_kmer_3bit.obj : encode_kmer_3bit.cpp encode_kmer_3bit.h

file.obj : file.cpp file.h file.h

asserts.obj : asserts.cpp asserts.h

clean :
	del fast_binary_search.exe *.obj
