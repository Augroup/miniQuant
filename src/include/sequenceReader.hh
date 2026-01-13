#ifndef SEQUENCE_READER_HH
#define SEQUENCE_READER_HH

#include <string>
#include <vector>

char getComplement(char nucleotide);
std::string get_reverse_complementary(const std::string& dnaSequence);
void reverse_complementary_for_all_sequences(std::vector<std::string>& reference_sequences,std::vector<std::string>& reverse_complement_sequences);
void readFASTA(const std::string& filename, std::vector<Transcript>& transcripts, std::vector<std::string>& sequenceData);
void getSequencesCharPtr(std::vector<std::string>& sequenceData, char*** sequences);
#endif