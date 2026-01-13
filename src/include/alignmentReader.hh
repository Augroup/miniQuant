#ifndef ALIGNMENT_READER_HH
#define ALIGNMENT_READER_HH

#include <string>
#include <vector>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include <common.hh>
#include <EM.hh>

// AlignmentReader class
class AlignmentReader {
public:
    // Method to read and process each alignment into Read and Hit structures directly
    void processLongReadsAlignments(std::vector<Transcript>& transcripts,
                                            std::vector<Read>& reads,
                                            std::vector<size_t>& thread_bucket_id,int chunk_id,
                                            ProgramOptions& opt);
    // Method to read and process each alignment into Read and Hit structures directly
    void processPairedShortReadsAlignments(
                                    std::vector<Transcript>& transcripts,
                                    std::vector<Read>& reads,
                                    std::vector<size_t>& thread_bucket_id,
                                    ProgramOptions& opt);
    void parallelProcessAlignments(
                                    std::vector<Transcript>& transcripts,
                                    std::vector<Read>& reads,
                                    ProgramOptions& opt,ReadType readType);
    
};

#endif // ALIGNMENTREADER_HH

