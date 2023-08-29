# CRISPR

This should compile and run under Linux, MacOS and Windows.  Its written using C++11.

Given a fasta file such as GCF_001433935.1_IRGSP-1.0_genomic.fna, its first necessary to extract the guides (queries) and off-site targets (documents).  The tool that does this is called extract_guides

```
Usage: ./extract_guides <fna_fileame> <output_file>
```

so

```
./extract_guides GCF_001433935.1_IRGSP-1.0_genomic.fna OryzaSativaGuides.txt
```


That file is formatted thus:

\<20mer\> \<strand\> \<sequence\> \<location\> \<sitetype\>

e.g.

```
AGGGTTTAGGGTTTAGGGTT  -  1 1001 A
```

where AGGGTTTAGGGTTTAGGGTT is the 20-mer, 1 is the chromosome, 1001 is the location on the chromosome, and A the that the N[AG]G at the end of the CRISPR site is NAG not NGG.


To do the search run fast_binary_search,

```
/Usage:./fast_binary_search [-b -h] [-f<filename>] [-o<filename>] [-t<threadcount>] [-debug]
       -? | -help print this help message
       -debug only search for the first 10,000 guides
       -b for binary search [default]
       -h for hamming distance
       -t<threads> the number of threads to use when searching [default = corecount (inc hyperthreads)]
       -f<filename> use <filename> as the genome
       -o<filename> use <filename> as the output file
```

where the default is the fast_binary_search and it will uses as many threads as you have cores/hyperthreads.

For example,

```
./fast_binary_search -f OryzaSativaGuides.txt -o OryzaSativaGuides.out.txt
```

Various statistics are reported including the time to search, and the "best" 20-mer.  The output fle reposts all the 20-mers that score sufficiently highly that that they might be useful (t <= 0.75 accumulated MIT score).


