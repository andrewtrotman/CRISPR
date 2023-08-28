/*
	EXTRACT_GUIDES.CPP
	------------------
*/
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <charconv>

#include "file.h"

uint8_t *chromosome = nullptr;

void reverse_complement_20mer(uint8_t *into, const uint8_t *from)
	{
	into += 19;
	const uint8_t *end = from + 20;
	while (from < end)
		{
		switch (*from)
			{
			case 'A':
				*into-- = 'T';
				break;
			case 'T':
				*into-- = 'A';
				break;
			case 'C':
				*into-- = 'G';
				break;
			case 'G':
				*into-- = 'C';
				break;
			default:
				*into-- = *from;
			}
		from++;
		}
	}

/*
	EXTRACT_20MERS()
	----------------
*/
void extract_20mers(JASS::file &into, const std::string sequence_number, const uint8_t *chromosome)
	{
	uint8_t flipped[40];
	char number[40];

	for (const uint8_t *base = chromosome + 21; *base != '\0'; base++)
		{
		/*
			According to
				J. Bradford, T. Chappell, D. Perrin (2022), Rapid Whole-Genome Identification of High Quality CRISPR Guide RNAs with the Crackling Method. The CRISPR Journal 5(3):410-421. http://doi.org/10.1089/crispr.2021.0102
			The sequence must not start with a T, and must end NAG or NGG (see section "Off-target scoring performance").
			That is, the regular expression is [ACG][ATCG]20[AG]G.  Were GG is a guid (query) and GG or AG is an off-target
			site (document).
		*/
		if (*(base + 1) == 'G' && (*base == 'A' || *base == 'G'))
			{
			if (*base - 21 == 'T')
				continue;
			/*
				If we have GG then its a guide (query), if we have AG or GG then its an off-target guide (document).
			*/
			bool failed = false;
			const uint8_t *stop_at = base - 20;
			for (const uint8_t *check = base; check > stop_at; check--)
				if (*check != 'A' && *check != 'C' && *check != 'G' && *check != 'T')
					{
					failed = true;
					break;
					}
			if (!failed)
				{
//				std::cout.write((const char *)base - 21, 20) <<  "  +  " << sequence_number << " " << base - chromosome - 20 << '\n';
				into.write((const char *)base - 21, 20);
				into.write("  +  ", 5);
				into.write(sequence_number.c_str(), sequence_number.size());
				into.write(" ", 1);

				auto [end_of_number, error_code] = std::to_chars(number, number + sizeof(number), base - chromosome - 20);
				*end_of_number++ = ' ';
				*end_of_number++ = *base;
				*end_of_number = '\n';

				into.write(number, end_of_number - number + 1);

				}
			}
		/*
			Include the reverse complement sequences.  That is, ones that the other DNA strand matches the regular expression [ACG][ATCG]20[AG]G.
			That is, this strand is C[TC][ACTG]20[TGC].

			NOTE: THIS CANNOT BE AN ELSE AS THE CC MIGHT OVERLAP A GG
		*/
		if (*(base - 20) == 'C' && (*(base - 19) == 'C' || *(base - 19) == 'T'))
			{
			if (*base == 'A')
				continue;

			/*
				If we have CC then its a guide (query), if we have TC or CC then its an off-target guide (document).
			*/
			bool failed = false;
			const uint8_t *stop_at = base + 3;
			for (const uint8_t *check = base - 17; check < stop_at; check++)
				if (*check != 'A' && *check != 'C' && *check != 'G' && *check != 'T')
					{
					failed = true;
					break;
					}
			if (!failed)
				{
//				std::cout << reverse_complement(std::string((const char *)base - 17, 20)) << "  -  " << sequence_number << " " << base - chromosome - 19 << '\n';
				reverse_complement_20mer(&flipped[0], base - 17);
				into.write(flipped, 20);
				into.write("  -  ", 5);
				into.write(sequence_number.c_str(), sequence_number.size());
				into.write(" ", 1);

				auto [end_of_number, error_code] = std::to_chars(number, number + sizeof(number), base - chromosome - 19);
				*end_of_number++ = ' ';
				*end_of_number++ = *(base - 19) == 'C' ? 'G' : 'A';				// complement of the start (C[TC] -> G[AG]) so output A for off-targets (documents) and G for guides (queries) or off-targets (documents).
				*end_of_number = '\n';
				into.write(number, end_of_number - number + 1);
				}
			}
		}
	}

/*
	PARSE_GENOME()
	--------------
*/
void parse_genome(const char *infilename, char *outfilename)
	{
	JASS::file::file_read_only genome;
	JASS::file outfile(outfilename, "wb");
//	outfile.setvbuf(128 * 1024 * 1024);


	size_t size = genome.open(infilename);
	if (size == 0)
		{
		std::cout << "Cannot open file " << infilename << "\n";
		exit(1);
		}

	const uint8_t *buffer;
	genome.read_entire_file(buffer);
	const uint8_t *end_of_file = buffer + size;

	const uint8_t *start = buffer;
	const uint8_t *end = start;
	uint64_t sequence_number = 0;
	std::string sequence_number_as_string;
	do
		{
		/*
			Search for the chromosome details
		*/
		start = end;
		while (*start != '>')
			start++;
		while (*end !='\n')
			end++;
		sequence_number++;
		std::cout.write((const char *)start, end - start) << '\n';
		sequence_number_as_string = std::to_string(sequence_number);
		/*
			Search for the end of the chromosome
		*/
		start = end;
		while (*end != '>' && end < end_of_file)
			end++;

		/*
			We can now copy it and remove the noise (while doing so)
		*/
		chromosome = (uint8_t *)realloc(chromosome, end - start + 1);
		const uint8_t *from = start;
		uint8_t *to = chromosome;
		while (from < end)
			{
			if (*from >= 'a' && *from <= 'z')
				*to++ = *from++ - 32;			// convert to uppercase
			else if (*from >= 'A' && *from <= 'Z')
				*to++ = *from++;
			else
				from++;
			}
		*to = '\0';		// NULL terminate the chromosome

		/*
			Now extract all the 20-mers that end GG
		*/
		if (to - chromosome > 20)
			extract_20mers(outfile, sequence_number_as_string, chromosome);
		}
	while (end < end_of_file);
	}

/*
	USAGE()
	-------
*/
int usage(char *exename)
	{
	std::cout << "Usage: " << exename << " <fna_fileame> <output_file>\n";
	return 1;
	}

/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	if (argc != 3)
		exit(usage(argv[0]));

	parse_genome(argv[1], argv[2]);

	return 0;
	}
