/*
	EXTRACT_GUIDES.CPP
	------------------
*/
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>

#include "file_read_only.h"

uint8_t *chromosome = nullptr;

std::string reverse_complement(std::string &&sequence)
	{
	std::string reverseSeq;
	for (char nucleotide : sequence)
		{
		switch (nucleotide)
			{
			case 'A':
				reverseSeq += 'T';
				break;
			case 'T':
				reverseSeq += 'A';
				break;
			case 'C':
				reverseSeq += 'G';
				break;
			case 'G':
				reverseSeq += 'C';
				break;
			default:
				reverseSeq += nucleotide;
			}
		}
	std::reverse(reverseSeq.begin(), reverseSeq.end());
	return reverseSeq;
	}

/*
	EXTRACT_20MERS()
	----------------
*/
void extract_20mers(uint64_t sequence_number, const uint8_t *chromosome)
	{
	for (const uint8_t *base = chromosome + 21; *base != '\0'; base++)
		{
		if (*base == 'G' && *(base + 1) == 'G')
			{
			bool failed = false;
			const uint8_t *stop_at = base - 20;
			const uint8_t *check;
			for (check = base; check > stop_at; check--)
				if (*check != 'A' && *check != 'C' && *check != 'G' && *check != 'T')
					{
					failed = true;
					break;
					}
			if (!failed)
				std::cout.write((const char *)base - 21, 20) <<  "  +  " << sequence_number << " " << base - chromosome - 20 << '\n';
			}
		if (*(base - 20) == 'C' && *(base - 19) == 'C')
			{
			bool failed = false;
			const uint8_t *stop_at = base;
			const uint8_t *check;
			for (check = base; check > stop_at; check++)
				if (*check != 'A' && *check != 'C' && *check != 'G' && *check != 'T')
					{
					failed = true;
					break;
					}
			if (!failed)
				std::cout << reverse_complement(std::string((const char *)base - 17, 20)) << "  -  " << sequence_number << " " << base - chromosome - 19 << '\n';

//				std::cout.write((const char *)base - 17, 20) << "  -  " << sequence_number << " " << base - chromosome - 19 << '\n';
			}
		}
	}

/*
	PARSE_GENOME()
	--------------
*/
void parse_genome(const char *filename)
	{
	file_read_only genome;

	size_t size = genome.open(filename);
	if (size == 0)
		{
		std::cout << "Cannot open file " << filename << "\n";
		exit(1);
		}

	const uint8_t *buffer;
	genome.read_entire_file(buffer);
	const uint8_t *end_of_file = buffer + size;

	const uint8_t *start = buffer;
	const uint8_t *end = start;
	uint64_t sequence_number = 0;
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
			extract_20mers(sequence_number, chromosome);

exit(0);
		}
	while (end < end_of_file);
	}

/*
	USAGE()
	-------
*/
int usage(char *exename)
	{
	std::cout << "Usage: " << exename << " <fna_fileame>\n";
	return 1;
	}

/*
	MAIN()
	------
*/
int main(int argc, char *argv[])
	{
	if (argc != 2)
		exit(usage(argv[0]));

	parse_genome(argv[1]);

	return 0;
	}
