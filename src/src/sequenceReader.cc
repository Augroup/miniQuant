
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <EM.hh>
#include <sequenceReader.hh>
char getComplement(char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A';
        case 'a': return 'T';
        case 't': return 'A';
        case 'c': return 'G';
        case 'g': return 'C';
        case 'u': return 'a';
        default: return 'N'; // N for any non-standard nucleotide
    }
}


void getSequencesCharPtr(std::vector<std::string>& sequenceData, char*** sequences){
    size_t num_transcripts = sequenceData.size();
    *sequences = new char*[num_transcripts];
    for (int i = 0; i < num_transcripts; ++i) {
         // Here we put the transcript ID into the first line of sequence. And parse it in the MemoryReader in align.cc
        std::string sequence_with_id = std::to_string(i)+'\n'+sequenceData[i];
        (*sequences)[i] = new char[sequence_with_id.length() + 1];
        // Copy std::string to char arrays
        std::strcpy((*sequences)[i], sequence_with_id.c_str());
    }
}
std::string get_reverse_complementary(const std::string& dnaSequence){
    std::string reverse_complement_sequence;
    reverse_complement_sequence.reserve(dnaSequence.length());    
    // Compute the complement
    for (auto it = dnaSequence.rbegin(); it != dnaSequence.rend(); ++it) {
        reverse_complement_sequence.push_back(getComplement(*it));
    }
    return reverse_complement_sequence;
}
void reverse_complementary_for_all_sequences(std::vector<std::string>& reference_sequences,std::vector<std::string> &reverse_complement_sequences){
    for (auto &dnaSequence : reference_sequences){
        reverse_complement_sequences.push_back(get_reverse_complementary(dnaSequence));
    }
}
void readFASTA(const std::string& filename, std::vector<Transcript>& transcripts, std::vector<std::string>& sequenceData) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string currentID;
    std::stringstream currentSeq;
    int counter = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') { // Header line
            if (!currentID.empty()) {
                // making the sequence to upper case
                std::string seq = currentSeq.str();
                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                transcripts.push_back(Transcript(counter,currentID,seq.length()));
                sequenceData.push_back(seq);
                currentSeq.str(""); // Clear the stringstream
                counter++;
            }
            currentID = line.substr(1); // Remove '>'
        } else { // Sequence line
            currentSeq << line;
        }
    }

    // Don't forget to add the last  sequence
    if (!currentID.empty()) {
        std::string seq = currentSeq.str();
        transcripts.push_back(Transcript(counter,currentID,seq.length()));
        sequenceData.push_back(seq);
    }

    file.close();
}